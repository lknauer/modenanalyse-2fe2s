@echo off
REM install.bat - Installs modenanalyse_2fe2s into the current Python environment.
REM v1.0.0: nisspec3 is no longer installed (see CHANGELOG.md, section 3.4).
REM         UMAP/HDBSCAN/SCSD are pulled in automatically as dependencies.

cd /d "%~dp0"

echo.
echo ===========================================================
echo   modenanalyse_2fe2s v1.0.0 - Installation
echo ===========================================================
echo.

set "WHL="
for %%f in (modenanalyse_2fe2s-*.whl) do set "WHL=%%f"

if not defined WHL (
    echo ERROR: modenanalyse_2fe2s-*.whl not found in %CD%
    pause
    exit /b 1
)
echo Found: %WHL%
echo.

REM Remove previous versions if present
echo Removing previous versions if present...
pip uninstall -y modenanalyse_2fe2s >nul 2>&1
pip uninstall -y nisspec3 >nul 2>&1
echo.

REM Fresh installation
echo Installing modenanalyse_2fe2s and dependencies...
pip install --upgrade "%WHL%"
if errorlevel 1 (
    echo.
    echo ERROR: Installation failed.
    pause
    exit /b 1
)

echo.
echo ===========================================================
echo   Installation successful.
echo ===========================================================
echo.

REM Version check
python -c "import modenanalyse_2fe2s; print(f'  modenanalyse_2fe2s v{modenanalyse_2fe2s.__version__}')"
python -c "import umap; print('  umap-learn OK')" 2>nul || echo   umap-learn missing
python -c "import hdbscan; print('  hdbscan OK')" 2>nul || echo   hdbscan missing
python -c "import scsdpy; print('  scsdpy OK')" 2>nul || echo   scsdpy not available (optional)
echo.
echo NOTE: For NIS spectra or Fe-pDOS, install nisspec3 separately
echo       (pip install nisspec3) and call it on the same log file.
echo.
echo First run:  modenanalyse-2fe2s --write-template run.toml
echo.
pause
