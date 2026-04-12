#!/bin/bash
# Install a macOS launchd agent that runs check_and_run.sh every 5 minutes.
# Run once: bash tools/install_daemon.sh
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
LABEL="com.deep-gw-cnv.autorun"
PLIST="$HOME/Library/LaunchAgents/$LABEL.plist"

cat > "$PLIST" <<EOF
<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE plist PUBLIC "-//Apple//DTD PLIST 1.0//EN"
    "http://www.apple.com/DTDs/PropertyList-1.0.dtd">
<plist version="1.0">
<dict>
    <key>Label</key>
    <string>$LABEL</string>
    <key>ProgramArguments</key>
    <array>
        <string>/bin/bash</string>
        <string>$REPO_ROOT/tools/check_and_run.sh</string>
    </array>
    <key>StartInterval</key>
    <integer>300</integer>
    <key>StandardOutPath</key>
    <string>$REPO_ROOT/tools/daemon.log</string>
    <key>StandardErrorPath</key>
    <string>$REPO_ROOT/tools/daemon.log</string>
    <key>RunAtLoad</key>
    <false/>
</dict>
</plist>
EOF

launchctl unload "$PLIST" 2>/dev/null || true
launchctl load "$PLIST"
echo "Daemon installed and running. Polls every 5 minutes."
echo "To uninstall: launchctl unload $PLIST && rm $PLIST"
