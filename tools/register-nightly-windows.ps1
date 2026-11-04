# Register the scheduled task that runs tools/nightly.sh at 02:00
# every day, as the current user. Safe to re-run: it replaces the task.
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
$script = (Resolve-Path (Join-Path $PSScriptRoot 'nightly.sh')).Path
# C:\a\b.sh -> /c/a/b.sh, the form Git Bash runs without path translation.
$unix = '/' + $script.Substring(0, 1).ToLowerInvariant() + ($script.Substring(2) -replace '\\', '/')

$action = New-ScheduledTaskAction -Execute $bash -Argument "--login -c '$unix'"
$trigger = New-ScheduledTaskTrigger -Daily -At '02:00'
$settings = New-ScheduledTaskSettingsSet -StartWhenAvailable -WakeToRun `
    -AllowStartIfOnBatteries -DontStopIfGoingOnBatteries `
    -ExecutionTimeLimit (New-TimeSpan -Hours 6) -MultipleInstances IgnoreNew
$user = "$env:USERDOMAIN\$env:USERNAME"
$task = @{ TaskName = $name; Action = $action; Trigger = $trigger; Settings = $settings; Force = $true
           Description = 'Fetch origin/main, build and publish the popkern Windows nightly (tools/nightly.sh).' }

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
Write-Host "logs:     $env:LOCALAPPDATA\popkern-nightly\logs"
