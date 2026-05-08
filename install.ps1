# install.ps1 - Installs modenanalyse_2fe2s into the current Python environment.
# UMAP/HDBSCAN/SCSD are installed automatically as dependencies.

Set-Location -Path $PSScriptRoot

Write-Host ""
Write-Host "==========================================================="
Write-Host "  modenanalyse_2fe2s - Installation" -ForegroundColor Cyan
Write-Host "==========================================================="
Write-Host ""

$mod = Get-ChildItem -Filter "modenanalyse_2fe2s-*.whl" | Select-Object -First 1
if (-not $mod) {
    Write-Host "ERROR: modenanalyse_2fe2s-*.whl not found in $PWD" -ForegroundColor Red
    pause
    exit 1
}
Write-Host "Found: $($mod.Name)"
Write-Host ""

# Remove old versions
Write-Host "Removing previous versions if present..."
pip uninstall -y modenanalyse_2fe2s 2>$null
Write-Host ""

# Fresh installation; dependencies will be pulled automatically from PyPI
Write-Host "Installing modenanalyse_2fe2s and dependencies..." -ForegroundColor Yellow
pip install --upgrade $mod.FullName
if ($LASTEXITCODE -ne 0) {
    Write-Host ""
    Write-Host "ERROR: Installation failed." -ForegroundColor Red
    pause
    exit 1
}

Write-Host ""
Write-Host "===========================================================" -ForegroundColor Green
Write-Host "  Installation successful." -ForegroundColor Green
Write-Host "===========================================================" -ForegroundColor Green
Write-Host ""

# Version check
python -c "import modenanalyse_2fe2s; print(f'  modenanalyse_2fe2s v{modenanalyse_2fe2s.__version__}')"
python -c "import umap"    2>$null; if ($LASTEXITCODE -eq 0) { Write-Host "  umap-learn       OK" } else { Write-Host "  umap-learn       MISSING" -ForegroundColor Yellow }
python -c "import hdbscan" 2>$null; if ($LASTEXITCODE -eq 0) { Write-Host "  hdbscan          OK" } else { Write-Host "  hdbscan          MISSING" -ForegroundColor Yellow }
python -c "import scsdpy"  2>$null; if ($LASTEXITCODE -eq 0) { Write-Host "  scsdpy           OK" } else { Write-Host "  scsdpy           MISSING" -ForegroundColor Yellow }
Write-Host ""
Write-Host "First run:"
Write-Host "  modenanalyse-2fe2s --write-template run.toml   (creates a template)"
Write-Host "  modenanalyse-2fe2s run.toml                    (runs the analysis)"
Write-Host ""
pause
