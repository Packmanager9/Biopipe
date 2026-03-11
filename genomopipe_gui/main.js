const { app, BrowserWindow, ipcMain, shell, dialog, Menu } = require('electron');
const path   = require('path');
const { spawn } = require('child_process');
const fs     = require('fs');
const os     = require('os');

let mainWindow;
let activePipeline = null;
let activeFeedback = null;

// ─── Persistent settings ──────────────────────────────────────────────────────
const SETTINGS_PATH = path.join(app.getPath('userData'), 'bioforge-settings.json');

const DEFAULT_SETTINGS = {
  condaProfile:  `${os.homedir()}/miniconda3/etc/profile.d/conda.sh`,
  condaBin:      `${os.homedir()}/miniconda3/bin/conda`,
  genemarkPath:  `${os.homedir()}/genemark-etp-full/gmetp_linux_64/bin`,
  scriptsDir:    __dirname,
  outputDir:     `${os.homedir()}/output`,
  jalviewCmd:    'jalview',
  snapgeneCmd:   'snapgene-viewer',
  pymolCmd:      'pymol',
  vmdCmd:        'vmd',
  chimeraCmd:    'chimerax',
};

function loadSettings() {
  try {
    if (fs.existsSync(SETTINGS_PATH))
      return { ...DEFAULT_SETTINGS, ...JSON.parse(fs.readFileSync(SETTINGS_PATH, 'utf8')) };
  } catch (_) {}
  return { ...DEFAULT_SETTINGS };
}

function saveSettings(s) {
  try { fs.writeFileSync(SETTINGS_PATH, JSON.stringify(s, null, 2)); return true; }
  catch (e) { return false; }
}

// ─── Window ───────────────────────────────────────────────────────────────────
function createWindow() {
  mainWindow = new BrowserWindow({
    width: 1480, height: 920, minWidth: 1100, minHeight: 700,
    backgroundColor: '#040e1c',
    titleBarStyle: process.platform === 'darwin' ? 'hiddenInset' : 'default',
    webPreferences: {
      preload: path.join(__dirname, 'preload.js'),
      nodeIntegration: false,
      contextIsolation: true,
    },
    show: false,
  });
  mainWindow.loadFile('index.html');
  mainWindow.once('ready-to-show', () => { mainWindow.show(); });
  mainWindow.on('closed', () => { mainWindow = null; });
}

app.on('ready', () => {
  Menu.setApplicationMenu(null);
  createWindow();
});
app.on('window-all-closed', () => { if (process.platform !== 'darwin') app.quit(); });
app.on('activate', () => { if (!mainWindow) createWindow(); });

// ─── Settings IPC ─────────────────────────────────────────────────────────────
ipcMain.handle('get-settings', ()      => loadSettings());
ipcMain.handle('save-settings', (_, s) => saveSettings(s));

// ─── File system IPC ──────────────────────────────────────────────────────────
ipcMain.handle('read-file', (_, p) => {
  try {
    if (!fs.existsSync(p)) return { ok: false, error: 'File not found' };
    const stat = fs.statSync(p);
    if (stat.size > 25 * 1024 * 1024) return { ok: false, error: 'File too large (>25 MB)' };
    return { ok: true, content: fs.readFileSync(p, 'utf8') };
  } catch (e) { return { ok: false, error: e.message }; }
});

ipcMain.handle('list-directory', (_, dir) => {
  try {
    if (!fs.existsSync(dir)) return { ok: false, error: 'Not found: ' + dir };
    const entries = fs.readdirSync(dir, { withFileTypes: true })
      .map(e => {
        const full = path.join(dir, e.name);
        let size = 0;
        try { if (e.isFile()) size = fs.statSync(full).size; } catch (_) {}
        return { name: e.name, path: full, isDir: e.isDirectory(), size };
      })
      .sort((a, b) => {
        if (a.isDir !== b.isDir) return a.isDir ? -1 : 1;
        return a.name.localeCompare(b.name);
      });
    return { ok: true, entries };
  } catch (e) { return { ok: false, error: e.message }; }
});

ipcMain.handle('path-exists', (_, p) => {
  try { return fs.existsSync(p); } catch (_) { return false; }
});

ipcMain.handle('get-file-stats', (_, p) => {
  try {
    const s = fs.statSync(p);
    return { ok: true, mtimeMs: s.mtimeMs, size: s.size };
  } catch (e) { return { ok: false, error: e.message }; }
});

ipcMain.handle('get-uptime-ms', () => {
  // os.uptime() returns seconds since last boot
  return os.uptime() * 1000;
});

ipcMain.handle('open-dir-dialog', async (_, defaultPath) => {
  let startPath = os.homedir();
  if (defaultPath && defaultPath.trim()) {
    try { if (fs.existsSync(defaultPath.trim())) startPath = defaultPath.trim(); } catch (_) {}
  }
  const r = await dialog.showOpenDialog(mainWindow, {
    defaultPath: startPath, properties: ['openDirectory'],
  });
  return r.canceled ? null : r.filePaths[0];
});

ipcMain.handle('open-file-dialog', async (_, opts) => {
  const r = await dialog.showOpenDialog(mainWindow, {
    defaultPath: opts?.defaultPath || os.homedir(),
    filters: opts?.filters || [{ name: 'All Files', extensions: ['*'] }],
    properties: ['openFile'],
  });
  return r.canceled ? null : r.filePaths[0];
});

// ─── Load structured results ──────────────────────────────────────────────────
ipcMain.handle('load-results', (_, inputDir) => {
  const resolveDir = (d) => {
    try {
      const link = path.join(d, 'latest');
      if (fs.existsSync(link) && fs.lstatSync(link).isSymbolicLink())
        return fs.realpathSync(link);
    } catch (_) {}
    return fs.existsSync(d) ? d : null;
  };

  const rd = resolveDir(inputDir);
  if (!rd) return { ok: false, error: 'Run directory not found: ' + inputDir };

  const read   = p => (fs.existsSync(p) ? fs.readFileSync(p, 'utf8') : null);
  const pdbs   = dir => !fs.existsSync(dir) ? [] :
    fs.readdirSync(dir).filter(f => f.endsWith('.pdb')).map(f => ({
      file: f, path: path.join(dir, f), content: fs.readFileSync(path.join(dir, f), 'utf8'),
    }));
  const fastas = dir => !fs.existsSync(dir) ? [] :
    fs.readdirSync(dir).filter(f => /\.(fasta|fa|faa)$/i.test(f)).map(f => ({
      file: f, path: path.join(dir, f), content: fs.readFileSync(path.join(dir, f), 'utf8'),
    }));
  const gbs    = dir => !fs.existsSync(dir) ? [] :
    fs.readdirSync(dir).filter(f => /\.(gb|gbk)$/i.test(f)).map(f => ({
      file: f, path: path.join(dir, f), content: fs.readFileSync(path.join(dir, f), 'utf8'),
    }));

  const sentinels = {};
  ['phase1','phase2','feedback1','feedback2','feedback3','feedback4','feedback5','feedback6']
    .forEach(s => { sentinels[s] = fs.existsSync(path.join(rd, `.genomopipe_${s}.done`)); });

  return {
    ok: true, runDir: rd, sentinels,
    proteins:    read(path.join(rd, 'proteins.faa')),
    designs:     pdbs(path.join(rd, 'designs')),
    structures:  pdbs(path.join(rd, 'colabfold_out')),
    splitSeqs:   fastas(path.join(rd, 'proteinmpnn_out', 'split_seqs')),
    plasmids:    gbs(path.join(rd, 'moclo_plasmids')),
    blastText:   read(path.join(rd, 'blast_results.txt')),
    pipelineLog: read(path.join(rd, 'logs', 'pipeline.log')),
    summaryMd:   read(path.join(rd, 'genomopipe_summary.md')),
    fb6Audit:    read(path.join(rd, 'feedback6_loop', 'feedback6_taxonomy_audit.txt')),
    fb4Report:   read(path.join(rd, 'feedback4_loop', 'fb4_report.tsv')),
    fb2Summary:  read(path.join(rd, 'feedback2_loop', 'fb2_summary.tsv')),
  };
});

// ─── Launch external tools ────────────────────────────────────────────────────
ipcMain.on('launch-tool', (event, { tool, filePath }) => {
  if (tool === 'reveal')       { shell.showItemInFolder(filePath); return; }
  if (tool === 'open-default') { shell.openPath(filePath); return; }

  const s = loadSettings();
  const CMD = { jalview: s.jalviewCmd, snapgene: s.snapgeneCmd,
                pymol: s.pymolCmd, vmd: s.vmdCmd, chimera: s.chimeraCmd };
  const cmd = CMD[tool];
  if (!cmd) return;

  const proc = spawn(cmd, filePath ? [filePath] : [], { detached: true, stdio: 'ignore' });
  proc.unref();
  proc.on('error', err => {
    if (err.code === 'ENOENT') shell.openPath(filePath);
    event.sender.send('tool-error', `${tool}: ${err.message}`);
  });
});

// ─── Run pipeline ─────────────────────────────────────────────────────────────
ipcMain.on('run-pipeline', (event, args) => {
  const s = loadSettings();
  const { organism, outputDir, isEukaryote, bam, autoRnaseq, force, configFile } = args;
  const scriptsDir = s.scriptsDir || '';
  const outDir     = outputDir || s.outputDir;

  // Resolve scripts — use full path if scriptsDir is set and file exists there,
  // otherwise fall back to bare name and let the shell resolve via PATH.
  const resolveScript = (name) => {
    if (scriptsDir) {
      const full = path.join(scriptsDir, name);
      if (fs.existsSync(full)) return { cmd: full, found: true };
    }
    return { cmd: name, found: false };
  };

  const orch   = resolveScript('genomopipe.py');
  const leg    = resolveScript('genome_to_design.sh');
  const useOrch = orch.found;

  const orchestrator = orch.cmd;
  const legacy       = leg.cmd;

  let fullCmd;
  if (useOrch) {
    const parts = [];
    if (configFile && fs.existsSync(configFile)) parts.push(`"${configFile}"`);
    if (organism)   parts.push(`--organism "${organism}"`);
    if (outDir)     parts.push(`--output_dir "${outDir}"`);
    if (isEukaryote) parts.push('--is_eukaryote true');
    if (bam)        parts.push(`--bam="${bam}"`);
    if (autoRnaseq) parts.push('--auto_rnaseq');
    if (force)      parts.push('--force');
    parts.push(`--scripts_dir "${scriptsDir}"`);
    parts.push(`--genemark_path "${s.genemarkPath}"`);
    fullCmd = `source "${s.condaProfile}" && conda activate braker_env && python "${orchestrator}" ${parts.join(' ')}`;
  } else {
    const parts = [`"${organism}"`, `"${outDir}"`];
    if (isEukaryote) parts.push('true');
    if (bam)         parts.push(`--bam="${bam}"`);
    if (autoRnaseq)  parts.push('--auto_rnaseq');
    if (force)       parts.push('--force');
    parts.push(`--GENEMARK_PATH="${s.genemarkPath}"`);
    fullCmd = `source "${s.condaProfile}" && conda activate braker_env && bash "${legacy}" ${parts.join(' ')}`;
  }

  event.sender.send('pipeline-log', `[BIOFORGE] Mode: ${useOrch ? 'orchestrator' : 'legacy'}\n`);
  event.sender.send('pipeline-log', `[CMD] ${fullCmd}\n`);

  const env = { ...process.env, PATH: `${process.env.PATH}:${s.genemarkPath}` };
  activePipeline = spawn('bash', ['-c', fullCmd], { env });
  activePipeline.stdout.on('data', d => event.sender.send('pipeline-log', d.toString()));
  activePipeline.stderr.on('data', d => event.sender.send('pipeline-log', d.toString()));
  activePipeline.on('close', code => {
    activePipeline = null;
    event.sender.send('pipeline-complete', { code, outputDir: outDir });
  });
  activePipeline.on('error', err =>
    event.sender.send('pipeline-log', `[ERROR] ${err.message}\n`));
});

ipcMain.on('kill-pipeline', () => {
  if (activePipeline) { activePipeline.kill('SIGTERM'); activePipeline = null; }
});

// ─── Run feedback loops ───────────────────────────────────────────────────────
ipcMain.on('run-feedback', (event, { loop, runDir, extraArgs }) => {
  const s = loadSettings();
  const scriptsDir = s.scriptsDir || '';

  const SCRIPTS = {
    fb1: ['feedback1_colabfold_to_rfdiffusion.sh', 'sh'],
    fb2: ['feedback2_plddt_mpnn_resample.py',       'py'],
    fb3: ['feedback3_blast_to_braker.sh',           'sh'],
    fb4: ['feedback4_domesticated_cds_revalidate.py','py'],
    fb5: ['feedback5_designed_proteins_to_annotation.sh','sh'],
    fb6: ['feedback6_blast_taxonomy_rerun.py',      'py'],
  };
  const entry = SCRIPTS[loop];
  if (!entry) { event.sender.send('feedback-log', `[ERROR] Unknown loop: ${loop}\n`); return; }

  const scriptPath = path.join(scriptsDir, entry[0]);
  if (!fs.existsSync(scriptPath)) {
    event.sender.send('feedback-log', `[ERROR] Script not found: ${scriptPath}\n`); return;
  }

  const argStr  = ['--run_dir', `"${runDir}"`, ...(extraArgs || [])].join(' ');
  const runner  = entry[1] === 'py' ? `python "${scriptPath}"` : `bash "${scriptPath}"`;
  const fullCmd = `source "${s.condaProfile}" && conda activate braker_env && ${runner} ${argStr}`;

  event.sender.send('feedback-log', `[BIOFORGE] ${loop}\n[CMD] ${fullCmd}\n`);
  if (activeFeedback) activeFeedback.kill('SIGTERM');
  activeFeedback = spawn('bash', ['-c', fullCmd]);
  activeFeedback.stdout.on('data', d => event.sender.send('feedback-log', d.toString()));
  activeFeedback.stderr.on('data', d => event.sender.send('feedback-log', d.toString()));
  activeFeedback.on('close', code => {
    activeFeedback = null;
    event.sender.send('feedback-complete', { loop, code });
  });
});

ipcMain.on('kill-feedback', () => {
  if (activeFeedback) { activeFeedback.kill('SIGTERM'); activeFeedback = null; }
});
