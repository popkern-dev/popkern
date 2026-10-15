# Prepare a Windows machine to build and publish popkern's Windows release
# archives (tools/dist-windows.sh). Safe to re-run: every step checks first
# and only installs what is missing.
#
#   powershell -ExecutionPolicy Bypass -File tools\setup-windows-build.ps1
#
# Installs, per user unless noted:
#   - Visual Studio 2022 Build Tools: the MSVC x64 and ARM64 compilers and
#     the Windows SDK (machine-wide; raises one UAC prompt)
#   - rustup, with the toolchain rust-toolchain.toml pins as the default and
#     the x86_64/aarch64 Windows and wasm32 targets
#   - wasm-bindgen-cli at the version Cargo.lock locks (tools/build-live.sh)
#   - Node.js LTS and wrangler (uploads to the popkern-releases R2 bucket)
# and checks that Git for Windows' bash is present (the tools are bash).
#
# It does not log wrangler in: run `wrangler login` once (a browser OAuth
# grant), or for unattended nightly builds set CLOUDFLARE_API_TOKEN (an R2
# edit token) and CLOUDFLARE_ACCOUNT_ID in the build user's environment.
$ErrorActionPreference = 'Stop'
$ProgressPreference = 'SilentlyContinue'
$repo = Split-Path -Parent $PSScriptRoot

function Step([string] $message) { Write-Host "== $message" }
function Have([string] $command) { [bool] (Get-Command $command -ErrorAction SilentlyContinue) }
# Pick up PATH changes installers just made, without a new terminal.
function Refresh-Path {
    $machine = [Environment]::GetEnvironmentVariable('Path', 'Machine')
    $user = [Environment]::GetEnvironmentVariable('Path', 'User')
    $env:Path = "$env:USERPROFILE\.cargo\bin;$machine;$user"
}
function Winget-Install([string] $id, [string[]] $extra = @()) {
    # One winget at a time: two concurrent installs corrupt its source
    # cache (0x80071130, "Fast Cache data not found").
    & winget install --id $id --exact --silent --accept-package-agreements --accept-source-agreements @extra
    if ($LASTEXITCODE -ne 0) { throw "winget install $id failed ($LASTEXITCODE)" }
}

if (-not (Have winget)) { throw 'winget is missing: install "App Installer" from the Microsoft Store' }

Step 'Visual Studio 2022 Build Tools (MSVC x64 + ARM64, Windows SDK)'
$vswhere = "${env:ProgramFiles(x86)}\Microsoft Visual Studio\Installer\vswhere.exe"
$components = @(
    'Microsoft.VisualStudio.Component.VC.Tools.x86.x64',
    'Microsoft.VisualStudio.Component.VC.Tools.ARM64',
    'Microsoft.VisualStudio.Component.Windows11SDK.22621'
)
$haveVs = (Test-Path $vswhere) -and (& $vswhere -products * -requires $components -property installationPath)
if ($haveVs) {
    Write-Host "   present: $haveVs"
} else {
    $add = ($components + 'Microsoft.VisualStudio.Workload.VCTools' | ForEach-Object { "--add $_" }) -join ' '
    Winget-Install 'Microsoft.VisualStudio.2022.BuildTools' @('--override', "--wait --passive --norestart $add")
}

Step 'Rust (rustup, pinned toolchain, release targets)'
Refresh-Path
if (-not (Have rustup)) {
    $dir = Join-Path $env:TEMP ("rustup-init-" + [Guid]::NewGuid().ToString('N'))
    New-Item -ItemType Directory -Path $dir | Out-Null
    try {
        Invoke-WebRequest 'https://static.rust-lang.org/rustup/dist/x86_64-pc-windows-msvc/rustup-init.exe' -OutFile "$dir\rustup-init.exe" -UseBasicParsing
        & "$dir\rustup-init.exe" -y --default-toolchain none --profile minimal
        if ($LASTEXITCODE -ne 0) { throw "rustup-init failed ($LASTEXITCODE)" }
    } finally {
        Remove-Item -Recurse -Force $dir -ErrorAction SilentlyContinue
    }
    Refresh-Path
}
$pinned = (Select-String -Path "$repo\rust-toolchain.toml" -Pattern '^channel = "(.+)"').Matches[0].Groups[1].Value
# The pin as the default too, so cargo also works outside the checkout
# (cargo install below) without fetching a second toolchain.
& rustup toolchain install $pinned --profile minimal --component rustfmt --component clippy
& rustup default $pinned
& rustup target add --toolchain $pinned x86_64-pc-windows-msvc aarch64-pc-windows-msvc wasm32-unknown-unknown
if ($LASTEXITCODE -ne 0) { throw 'rustup target add failed' }

Step 'wasm-bindgen-cli (the version Cargo.lock locks)'
$lock = Get-Content "$repo\Cargo.lock"
$at = [Array]::IndexOf($lock, 'name = "wasm-bindgen"')
if ($at -lt 0) { throw 'cannot find wasm-bindgen in Cargo.lock' }
$locked = ($lock[$at + 1] -replace '^version = "(.+)"$', '$1')
$current = if (Have wasm-bindgen) { ((& wasm-bindgen --version) -split ' ')[1] } else { $null }
if ($current -eq $locked) {
    Write-Host "   present: wasm-bindgen $current"
} else {
    Push-Location $env:USERPROFILE
    try {
        & cargo install wasm-bindgen-cli --version $locked --locked --force
        if ($LASTEXITCODE -ne 0) { throw "cargo install wasm-bindgen-cli $locked failed" }
    } finally { Pop-Location }
}

Step 'Node.js LTS and wrangler'
Refresh-Path
if (-not (Have node)) { Winget-Install 'OpenJS.NodeJS.LTS'; Refresh-Path }
if (-not (Have wrangler)) {
    & npm install --global wrangler
    if ($LASTEXITCODE -ne 0) { throw 'npm install --global wrangler failed' }
    Refresh-Path
}
Write-Host "   wrangler $(& wrangler --version)"

Step 'Git for Windows bash'
$bash = Join-Path $env:ProgramFiles 'Git\bin\bash.exe'
if (-not (Test-Path $bash)) { Winget-Install 'Git.Git' }
Write-Host "   $bash"

Write-Host @'

Ready. Next, once:  wrangler login
Then, from Git Bash in the checkout:
    tools/dist-windows.sh                       # build + verify the archives
    tools/dist-windows.sh --upload [--latest]   # publish to dl.popkern.dev
'@
