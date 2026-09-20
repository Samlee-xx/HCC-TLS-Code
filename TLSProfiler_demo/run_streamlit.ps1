param(
    [string]$PythonExe = "python"
)

$ErrorActionPreference = "Stop"
$ProjectRoot = Split-Path -Parent $MyInvocation.MyCommand.Path
Set-Location $ProjectRoot

if (-not $env:DEEPCELL_ACCESS_TOKEN) {
    throw "DEEPCELL_ACCESS_TOKEN is not set in this PowerShell session."
}

& $PythonExe -m streamlit run (Join-Path $ProjectRoot "app.py")
exit $LASTEXITCODE
