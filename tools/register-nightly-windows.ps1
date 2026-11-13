# Register the scheduled task that runs tools/nightly.sh at 02:00
# every day, as the current user. Safe to re-run: it replaces the task.
#
# The task runs an installed copy of tools/nightly-bootstrap.sh
# (%LOCALAPPDATA%\popkern-nightly\bootstrap.sh), which runs origin/main's
# tools/nightly.sh from a runner worktree: a night never depends on what
# this checkout has checked out. Re-run this script when the bootstrap
# changes on main (a night's log says so); the checkout must stay where it
# is, since the runner and build worktrees share its repository.
#
#   powershell -ExecutionPolicy Bypass -File tools\register-nightly-windows.ps1
#   powershell -ExecutionPolicy Bypass -File tools\register-nightly-windows.ps1 -Unregister
#
# The task survives restarts (Task Scheduler keeps it), wakes the machine
# from sleep to run, and runs at the next boot or wake when 02:00 passed
# while the machine was off. Run from an elevated PowerShell, it runs
# whether or not you are logged on (S4U: no password is stored). Without
# elevation it can only run while you are logged on, and says so.
#
# Prerequisites: tools/setup-windows-build.ps1, then `wrangler login` once.
param([switch] $Unregister)
$ErrorActionPreference = 'Stop'
$name = 'popkern nightly Windows build'

if ($Unregister) {
    Unregister-ScheduledTask -TaskName $name -Confirm:$false
    Write-Host "removed task '$name'"
    return
}

$bash = Join-Path $env:ProgramFiles 'Git\bin\bash.exe'
if (-not (Test-Path $bash)) { throw "missing $bash (run tools\setup-windows-build.ps1)" }
# C:\a\b.sh -> /c/a/b.sh, the form Git Bash runs without path translation.
function Unix([string] $path) { '/' + $path.Substring(0, 1).ToLowerInvariant() + ($path.Substring(2) -replace '\\', '/') }
$repo = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
# Where tools/nightly.sh keeps its worktrees and logs (its own default).
$nightlyHome = if ($env:POPKERN_NIGHTLY_HOME) { $env:POPKERN_NIGHTLY_HOME } else { Join-Path $env:LOCALAPPDATA 'popkern-nightly' }
New-Item -ItemType Directory -Force -Path (Join-Path $nightlyHome 'logs') | Out-Null
$bootstrap = Join-Path $nightlyHome 'bootstrap.sh'
Copy-Item -LiteralPath (Join-Path $PSScriptRoot 'nightly-bootstrap.sh') -Destination $bootstrap -Force

$action = New-ScheduledTaskAction -Execute $bash -Argument "--login -c `"'$(Unix $bootstrap)' '$(Unix $repo)'`""
$trigger = New-ScheduledTaskTrigger -Daily -At '02:00'
$settings = New-ScheduledTaskSettingsSet -StartWhenAvailable -WakeToRun `
    -AllowStartIfOnBatteries -DontStopIfGoingOnBatteries `
    -ExecutionTimeLimit (New-TimeSpan -Hours 6) -MultipleInstances IgnoreNew
$user = "$env:USERDOMAIN\$env:USERNAME"
$task = @{ TaskName = $name; Action = $action; Trigger = $trigger; Settings = $settings; Force = $true
           Description = 'Build and publish the popkern Windows nightly with origin/main''s tools/nightly.sh (via bootstrap.sh).' }

try {
    $principal = New-ScheduledTaskPrincipal -UserId $user -LogonType S4U -RunLevel Limited
    Register-ScheduledTask @task -Principal $principal | Out-Null
    Write-Host "registered '$name': daily 02:00 as $user, whether or not logged on"
} catch {
    $principal = New-ScheduledTaskPrincipal -UserId $user -LogonType Interactive -RunLevel Limited
    Register-ScheduledTask @task -Principal $principal | Out-Null
    Write-Host "registered '$name': daily 02:00 as $user, ONLY while logged on"
    Write-Host '(re-run from an elevated PowerShell to run it while logged off)'
}
Write-Host "run now:  Start-ScheduledTask -TaskName '$name'"
Write-Host "logs:     $nightlyHome\logs"
