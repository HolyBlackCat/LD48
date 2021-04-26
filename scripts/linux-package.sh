set -euxo pipefail
rm -rf "LD48_earthwyrm_v${VER}_linux"
rm -rf "LD48_earthwyrm_v${VER}_linux.zip"
cp -af bin "LD48_earthwyrm_v${VER}_linux"
cd "LD48_earthwyrm_v${VER}_linux"
rm -f progress.txt
find . -name '_*' -prune -print -exec rm -rf "{}" \;
cd ..
tar -czf "LD48_earthwyrm_v${VER}_linux.tar.gz" "LD48_earthwyrm_v${VER}_linux"
rm -rf "LD48_earthwyrm_v${VER}_linux"
