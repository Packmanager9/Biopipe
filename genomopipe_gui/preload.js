const { contextBridge, ipcRenderer } = require('electron');

contextBridge.exposeInMainWorld('api', {
  // Settings
  getSettings:  ()  => ipcRenderer.invoke('get-settings'),
  saveSettings: (s) => ipcRenderer.invoke('save-settings', s),

  // File system
  readFile:       (p)    => ipcRenderer.invoke('read-file', p),
  listDirectory:  (p)    => ipcRenderer.invoke('list-directory', p),
  pathExists:     (p)    => ipcRenderer.invoke('path-exists', p),
  getFileStats:   (p)    => ipcRenderer.invoke('get-file-stats', p),
  getUptimeMs:    ()     => ipcRenderer.invoke('get-uptime-ms'),
  openDirDialog:  (p)    => ipcRenderer.invoke('open-dir-dialog', p),
  openFileDialog: (opts) => ipcRenderer.invoke('open-file-dialog', opts),
  loadResults:    (dir)  => ipcRenderer.invoke('load-results', dir),

  // External tools
  launchTool:  (tool, filePath) => ipcRenderer.send('launch-tool', { tool, filePath }),
  onToolError: (cb) => ipcRenderer.on('tool-error', (_, msg) => cb(msg)),

  resumeDetect: (dir) => ipcRenderer.invoke('resume-detect', dir),

  // Pipeline
  runPipeline:    (args) => ipcRenderer.send('run-pipeline', args),
  killPipeline:   ()     => ipcRenderer.send('kill-pipeline'),
  onPipelineLog:  (cb)   => ipcRenderer.on('pipeline-log',     (_, d) => cb(d)),
  onPipelineDone: (cb)   => ipcRenderer.on('pipeline-complete', (_, d) => cb(d)),

  // Feedback loops
  runFeedback:    (args) => ipcRenderer.send('run-feedback', args),
  killFeedback:   ()     => ipcRenderer.send('kill-feedback'),
  onFeedbackLog:  (cb)   => ipcRenderer.on('feedback-log',     (_, d) => cb(d)),
  onFeedbackDone: (cb)   => ipcRenderer.on('feedback-complete', (_, d) => cb(d)),
});
