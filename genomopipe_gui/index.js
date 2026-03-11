// main.js - Electron Main Process
const { app, BrowserWindow, ipcMain } = require('electron');
const path = require('path');
const { spawn } = require('child_process');
const fs = require('fs');

let mainWindow;

function createWindow() {
  mainWindow = new BrowserWindow({
    width: 1280,
    height: 720,
    webPreferences: {
      preload: path.join(__dirname, 'preload.js'),
      nodeIntegration: false,
      contextIsolation: true,
    }
  });

  mainWindow.loadFile('index.html');

  mainWindow.on('closed', function () {
    mainWindow = null;
  });
}

app.on('ready', createWindow);

app.on('window-all-closed', function () {
  if (process.platform !== 'darwin') app.quit();
});

app.on('activate', function () {
  if (mainWindow === null) createWindow();
});

// IPC to run the pipeline
ipcMain.on('run-pipeline', (event, args) => {
  const { organism, outputDir, isEukaryote, bam, autoRnaseq, force } = args;
  const scriptPath = path.join(__dirname, 'genome_to_design.sh'); // Assume script is in app dir

  let cmdArgs = [organism, outputDir];
  if (isEukaryote) cmdArgs.push('true');
  if (bam) cmdArgs.push(`--bam=${bam}`);
  if (autoRnaseq) cmdArgs.push('--auto_rnaseq');
  if (force) cmdArgs.push('--force');

  const process = spawn('bash', [scriptPath, ...cmdArgs]);

  process.stdout.on('data', (data) => {
    event.sender.send('pipeline-log', data.toString());
  });

  process.stderr.on('data', (data) => {
    event.sender.send('pipeline-log', data.toString());
  });

  process.on('close', (code) => {
    event.sender.send('pipeline-complete', { code, outputDir });
  });
});

// IPC to load results for visualization
ipcMain.on('load-results', (event, outputDir) => {
  // Parse key files: proteins.faa (FASTA), designs/*.pdb, colabfold_out/*.pdb, blast_results.txt
  const results = {
    proteins: parseFasta(path.join(outputDir, 'proteins.faa')),
    designs: parsePdbs(path.join(outputDir, 'designs')),
    structures: parsePdbs(path.join(outputDir, 'colabfold_out')),
    blast: fs.readFileSync(path.join(outputDir, 'blast_results.txt'), 'utf8') || ''
  };
  event.sender.send('results-loaded', results);
});

// Simple FASTA parser
function parseFasta(filePath) {
  if (!fs.existsSync(filePath)) return [];
  const content = fs.readFileSync(filePath, 'utf8');
  const sequences = [];
  let current = { id: '', seq: '' };
  content.split('\n').forEach(line => {
    if (line.startsWith('>')) {
      if (current.id) sequences.push(current);
      current = { id: line.slice(1), seq: '' };
    } else {
      current.seq += line.trim();
    }
  });
  if (current.id) sequences.push(current);
  return sequences;
}

// Simple PDB parser: Extract ATOM coordinates for visualization
function parsePdbs(dirPath) {
  if (!fs.existsSync(dirPath)) return [];
  const files = fs.readdirSync(dirPath).filter(f => f.endsWith('.pdb'));
  return files.map(file => {
    const content = fs.readFileSync(path.join(dirPath, file), 'utf8');
    const atoms = content.split('\n')
      .filter(line => line.startsWith('ATOM'))
      .map(line => {
        const x = parseFloat(line.slice(30, 38).trim());
        const y = parseFloat(line.slice(38, 46).trim());
        const z = parseFloat(line.slice(46, 54).trim());
        return { x, y, z };
      });
    return { file, atoms };
  });
}