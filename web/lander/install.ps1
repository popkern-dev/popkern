# Install popkern on Windows: the `popkern` command and popkern live's page.
#
#   irm https://popkern.dev/install.ps1 | iex
#
# What it does, in order:
#   1. picks the build for this machine (x86_64, or ARM64 when the release
#      has one; x86_64 runs on ARM64 Windows 11 under emulation otherwise).
#      The binary links the C runtime statically, so it needs no Visual C++
#      redistributable;
#   2. downloads it and the release's SHA256SUMS, and refuses to install
#      anything whose checksum does not match;
#   3. unpacks it to $POPKERN_HOME\<version> (default
#      %LOCALAPPDATA%\popkern), points the $POPKERN_HOME\current junction at
#      it, and adds $POPKERN_HOME\current\bin to your user PATH.
# Nothing outside that folder and your user PATH is touched, and nothing
# needs administrator rights.
#
# Settings (environment variables, set before running):
#   POPKERN_VERSION         a release to install, e.g. 0.1.0 (default: latest)
#   POPKERN_HOME            where releases unpack (default %LOCALAPPDATA%\popkern)
#   POPKERN_DOWNLOAD_BASE   where releases are fetched from
#                           (default https://dl.popkern.dev)
#   POPKERN_NO_MODIFY_PATH  set to 1 to leave your PATH alone
#
# Uninstall: delete %LOCALAPPDATA%\popkern and remove
# %LOCALAPPDATA%\popkern\current\bin from your user PATH
# (Settings > System > About > Advanced system settings > Environment Variables).
#
# Runs on Windows PowerShell 5.1 and PowerShell 7. The script is ASCII
# only: `irm` decodes a response without a charset as Latin-1. Everything
# lives in one function, so the script runs only once all of it has
# arrived, and it never calls `exit`, which would close the window that
# ran `iex`.

function Install-Popkern {
    $ErrorActionPreference = 'Stop'
    # Windows PowerShell 5.1 redraws its progress bar for every chunk, which
    # makes a download many times slower.
    $ProgressPreference = 'SilentlyContinue'
    Set-StrictMode -Version 3.0

    function Say([string] $message) { Write-Host "popkern: $message" }

    $base = if ($env:POPKERN_DOWNLOAD_BASE) { $env:POPKERN_DOWNLOAD_BASE } else { 'https://dl.popkern.dev' }
    $base = $base.TrimEnd('/')
    $popkernHome = if ($env:POPKERN_HOME) { $env:POPKERN_HOME } else { Join-Path $env:LOCALAPPDATA 'popkern' }

    if ([Environment]::OSVersion.Platform -ne [PlatformID]::Win32NT) {
        throw "popkern: install.ps1 is for Windows. On macOS and Linux run:  curl -fsSL https://popkern.dev/install.sh | sh"
    }
    if (-not [Environment]::Is64BitOperatingSystem) {
        throw "popkern: no build for 32-bit Windows. Try popkern live in your browser: https://popkern.dev/live/"
    }

    # The OS's own architecture, not this process's: a 32-bit or emulated
    # x64 PowerShell on ARM64 reports its own.
    $arch = $null
    try {
        $arch = [string] [System.Runtime.InteropServices.RuntimeInformation]::OSArchitecture
    } catch {
        $arch = if ($env:PROCESSOR_ARCHITEW6432) { $env:PROCESSOR_ARCHITEW6432 } else { $env:PROCESSOR_ARCHITECTURE }
    }
    # Best first. @() keeps a one-element answer an array (a bare switch
    # unrolls it to a string, and [0] would be its first character).
    $targets = @(switch -Regex ($arch) {
        '^(X64|AMD64)$' { 'x86_64-pc-windows-msvc' }
        '^Arm64$' { 'aarch64-pc-windows-msvc'; 'x86_64-pc-windows-msvc' }
        default { throw "popkern: no Windows build for $arch yet" }
    })

    # Windows PowerShell 5.1 on older .NET still offers TLS 1.0 first.
    try {
        [Net.ServicePointManager]::SecurityProtocol = [Net.ServicePointManager]::SecurityProtocol -bor [Net.SecurityProtocolType]::Tls12
    } catch {}

    $tmp = Join-Path ([IO.Path]::GetTempPath()) ("popkern-install-" + [Guid]::NewGuid().ToString('N'))
    New-Item -ItemType Directory -Path $tmp | Out-Null
    try {
        function Fetch([string] $url, [string] $file) {
            Invoke-WebRequest -Uri $url -OutFile $file -UseBasicParsing
        }

        $version = $env:POPKERN_VERSION
        if (-not $version) {
            try { Fetch "$base/latest" "$tmp\latest" } catch { throw "popkern: cannot reach $base/latest ($($_.Exception.Message))" }
            $version = ([IO.File]::ReadAllText("$tmp\latest")).Trim()
        }
        $version = $version -replace '^v', ''
        if ($version -notmatch '^[0-9A-Za-z.+-]+$') {
            throw "popkern: unexpected version '$version' from $base/latest"
        }

        try { Fetch "$base/v$version/SHA256SUMS" "$tmp\SHA256SUMS" } catch { throw "popkern: release $version has no SHA256SUMS ($base/v$version/)" }
        $sums = @{}
        foreach ($line in [IO.File]::ReadAllLines("$tmp\SHA256SUMS")) {
            if ($line -match '^([0-9a-fA-F]{64}) [ *](.+)$') { $sums[$Matches[2].Trim()] = $Matches[1].ToLowerInvariant() }
        }

        # The first build this machine can run that the release lists.
        $target = $null
        foreach ($candidate in $targets) {
            if ($sums.ContainsKey("popkern-$version-$candidate.zip")) { $target = $candidate; break }
        }
        if (-not $target) { throw "popkern: release $version has no Windows build for $arch ($base/v$version/)" }
        if ($target -ne $targets[0]) { Say "no ARM64 build in $version; installing the x86_64 build (it runs under emulation)" }

        $name = "popkern-$version-$target"
        $archive = "$name.zip"
        Say "downloading popkern $version ($target)"
        try { Fetch "$base/v$version/$archive" "$tmp\$archive" } catch { throw "popkern: no $archive in release $version ($base/v$version/)" }

        $want = $sums[$archive]
        $got = (Get-FileHash -Algorithm SHA256 -LiteralPath "$tmp\$archive").Hash.ToLowerInvariant()
        if ($got -ne $want) { throw "popkern: checksum mismatch for $archive (expected $want, got $got)" }

        Add-Type -AssemblyName System.IO.Compression.FileSystem
        [IO.Compression.ZipFile]::ExtractToDirectory("$tmp\$archive", "$tmp\unpack")
        $unpacked = "$tmp\unpack\$name"
        if (-not (Test-Path -LiteralPath "$unpacked\bin\popkern.exe" -PathType Leaf)) { throw "popkern: $archive has no bin\popkern.exe" }

        New-Item -ItemType Directory -Force -Path $popkernHome | Out-Null
        $dest = Join-Path $popkernHome $version
        $staged = "$dest.new"
        # Each replaced copy gets its own name: one still held open by a
        # running popkern cannot be deleted yet, and must not block the next
        # install. Earlier ones are swept here once nothing holds them.
        $retired = "$dest.old-" + [Guid]::NewGuid().ToString('N').Substring(0, 8)
        $leftovers = @(Get-ChildItem -LiteralPath $popkernHome -Directory -Filter "$version.old-*" -ErrorAction SilentlyContinue | ForEach-Object { $_.FullName })
        foreach ($leftover in @($staged) + $leftovers) {
            if (Test-Path -LiteralPath $leftover) { Remove-Item -LiteralPath $leftover -Recurse -Force -ErrorAction SilentlyContinue }
        }
        Move-Item -LiteralPath $unpacked -Destination $staged
        if (Test-Path -LiteralPath $dest) {
            # Reinstalling a version: move the old copy aside in one rename,
            # which either happens or changes nothing. (Deleting it in place
            # would remove files until it hit the open popkern.exe of a running
            # `popkern serve`, leaving a half-deleted release behind.) Windows
            # does allow renaming the folder of a running exe, so a running
            # `popkern serve` keeps going from the retired copy.
            try { [IO.Directory]::Move($dest, $retired) }
            catch {
                Remove-Item -LiteralPath $staged -Recurse -Force -ErrorAction SilentlyContinue
                throw "popkern: cannot replace $dest; close any running popkern (popkern serve) and try again"
            }
        }
        Move-Item -LiteralPath $staged -Destination $dest
        if (Test-Path -LiteralPath $retired) { Remove-Item -LiteralPath $retired -Recurse -Force -ErrorAction SilentlyContinue }

        # `current` is a directory junction (no administrator rights or
        # developer mode needed, unlike a symlink). Delete only the junction
        # itself: a recursive Remove-Item follows it on Windows PowerShell 5.1
        # and empties the release it points at.
        $current = Join-Path $popkernHome 'current'
        if (Test-Path -LiteralPath $current) { [IO.Directory]::Delete($current, $false) }
        New-Item -ItemType Junction -Path $current -Target $dest | Out-Null

        $binDir = Join-Path $current 'bin'
        $exe = Join-Path $binDir 'popkern.exe'
        $installed = $null
        try {
            $output = & $exe --version
            if ($LASTEXITCODE -eq 0) { $installed = $output }
        } catch {}
        if (-not $installed) {
            throw "popkern: installed, but $exe does not run on this machine"
        }
        Say "installed $installed to $dest"

        $onPath = @($env:Path -split ';' | Where-Object { $_ -and ($_.TrimEnd('\') -ieq $binDir) }).Count -gt 0
        if ($env:POPKERN_NO_MODIFY_PATH -eq '1') {
            if (-not $onPath) { Say "$binDir is not on your PATH (POPKERN_NO_MODIFY_PATH is set, so it was left alone)" }
        } else {
            # Read and write the raw registry value: Environment.GetEnvironmentVariable
            # expands %VARS% and SetEnvironmentVariable writes REG_SZ, which
            # together would freeze every %VAR% in the user's PATH.
            $key = [Microsoft.Win32.Registry]::CurrentUser.OpenSubKey('Environment', $true)
            try {
                $userPath = [string] $key.GetValue('Path', '', [Microsoft.Win32.RegistryValueOptions]::DoNotExpandEnvironmentNames)
                $entries = @($userPath -split ';' | Where-Object { $_ })
                if (-not ($entries | Where-Object { $_.TrimEnd('\') -ieq $binDir })) {
                    $key.SetValue('Path', (($entries + $binDir) -join ';'), [Microsoft.Win32.RegistryValueKind]::ExpandString)
                    # Tell Explorer (and new terminals it starts) that the
                    # environment changed: setting any user variable broadcasts
                    # WM_SETTINGCHANGE.
                    [Environment]::SetEnvironmentVariable('POPKERN_INSTALL_PATH_REFRESH', '1', 'User')
                    [Environment]::SetEnvironmentVariable('POPKERN_INSTALL_PATH_REFRESH', $null, 'User')
                    Say "added $binDir to your user PATH"
                }
            } finally {
                $key.Close()
            }
            if (-not $onPath) { $env:Path = "$env:Path;$binDir" }
        }

        Write-Host @'

Make your first part:

    popkern new my-parts/bracket
    popkern serve my-parts

then open http://127.0.0.1:7690 in your browser.
(Other terminals that were already open need a restart to find popkern.)
'@
    } finally {
        Remove-Item -LiteralPath $tmp -Recurse -Force -ErrorAction SilentlyContinue
    }
}

Install-Popkern
