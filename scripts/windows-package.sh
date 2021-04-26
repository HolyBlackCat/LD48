set -euxo pipefail
rm -rf "LD48_earthwyrm_v${VER}_windows"
rm -rf "LD48_earthwyrm_v${VER}_windows.zip"
cp -af win_bin "LD48_earthwyrm_v${VER}_windows"
cd "LD48_earthwyrm_v${VER}_windows"
rm -f progress.txt
find . -name '_*' -prune -print -exec rm -rf "{}" \;
cd ..
zip -rq "LD48_earthwyrm_v${VER}_windows.zip" "LD48_earthwyrm_v${VER}_windows"
rm -rf "LD48_earthwyrm_v${VER}_windows"
