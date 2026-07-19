param(
    [switch]$Full
)

$ErrorActionPreference = "Stop"
$repo = Resolve-Path (Join-Path $PSScriptRoot "..\..")
$stamp = Get-Date -Format "yyyyMMdd_HHmmss"
$runsRoot = Join-Path $repo "test_e2e_2026\99_runs"
New-Item -ItemType Directory -Path $runsRoot -Force -ErrorAction Stop | Out-Null
$runDir = Join-Path $runsRoot "final_release_$stamp"
New-Item -ItemType Directory -Path $runDir -Force -ErrorAction Stop | Out-Null

$summary = @()

function Run-Step {
    param(
        [string]$Name,
        [string]$Command
    )

    $safeName = $Name -replace '[^A-Za-z0-9_.-]', '_'
    $stdout = Join-Path $runDir "$safeName.stdout.log"
    $stderr = Join-Path $runDir "$safeName.stderr.log"
    $codeFile = Join-Path $runDir "$safeName.exit_code.txt"
    $started = Get-Date
    $code = 0

    Push-Location $repo
    try {
        & powershell -NoProfile -ExecutionPolicy Bypass -Command $Command > $stdout 2> $stderr
        $code = $LASTEXITCODE
        if ($null -eq $code) { $code = 0 }
    }
    catch {
        $_ | Out-File -FilePath $stderr -Append -Encoding UTF8
        $code = 1
    }
    finally {
        Pop-Location
    }

    Set-Content -Path $codeFile -Value $code -Encoding UTF8
    $elapsed = [Math]::Round(((Get-Date) - $started).TotalSeconds, 2)
    $status = if ($code -eq 0) { "PASS" } else { "FAIL" }
    $script:summary += [pscustomobject]@{
        name = $Name
        status = $status
        exit_code = $code
        seconds = $elapsed
        stdout = Split-Path $stdout -Leaf
        stderr = Split-Path $stderr -Leaf
    }
    Write-Host "$status $Name ($elapsed s)"
}

Run-Step "git_status_clean" "`$s = git status --porcelain; if (`$s) { `$s; exit 1 }"
Run-Step "python_compile" "python -m py_compile allenricher/cli.py allenricher/api/server.py allenricher/report/generator.py allenricher/ai/interpreter.py allenricher/core/enrichment.py allenricher/core/gsva.py allenricher/database/species_registry.py"
Run-Step "pytest_release_focus" "python -m pytest tests/test_cli.py tests/test_api_server.py tests/test_ai_interpreter.py tests/test_ai_integration.py tests/test_ai_structured.py tests/test_report_integration.py tests/test_web_database_defaults.py tests/test_cli_species_registry_parity.py tests/test_database_manager_correctness.py tests/test_database_manager_tf.py tests/test_method_correctness.py -q --override-ini addopts='' -p no:cacheprovider"
Run-Step "cli_help_smoke" "python -m allenricher --help"
Run-Step "cli_species_smoke" "python -m allenricher query-species --taxid 9606"

if ($Full) {
    Run-Step "pytest_full" "python -m pytest -q --override-ini addopts='' -p no:cacheprovider"
    Run-Step "real_world_sci_matrix" "python test_e2e_2026\18_real_world_sci\run_real_world_sci.py --mode local --keep-going"
}

$summaryPath = Join-Path $runDir "FINAL_RELEASE_SUMMARY.tsv"
$summary | Export-Csv -Path $summaryPath -NoTypeInformation -Delimiter "`t" -Encoding UTF8

$failed = @($summary | Where-Object { $_.status -ne "PASS" })
$decision = if ($failed.Count -eq 0) { "GO" } else { "NO_GO" }
Set-Content -Path (Join-Path $runDir "FINAL_RELEASE_DECISION.txt") -Value $decision -Encoding UTF8

Write-Host ""
Write-Host "Release decision: $decision"
Write-Host "Evidence directory: $runDir"

if ($failed.Count -gt 0) {
    Write-Host "Failed steps:"
    $failed | ForEach-Object { Write-Host "- $($_.name): exit $($_.exit_code)" }
    exit 1
}
