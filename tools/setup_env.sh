#!/usr/bin/env bash
# Provisions both stacks into project-local directories: .venv/ and .rlib/. See REPRODUCIBILITY.md, Environment layout and overrides, for the rationale.
#
# Usage: tools/setup_env.sh [--python-only|--r-only]
#
# Env:
#   IGNITE_CXX_SDK_FALLBACK=1   Builds against the SDK's libc++ when a probe confirms it compiles; opt-in fallback for a broken default toolchain.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

VENV="$ROOT/.venv"
RLIB="$ROOT/.rlib"
DO_PY=1; DO_R=1
case "${1:-}" in
  --python-only) DO_R=0 ;;
  --r-only)      DO_PY=0 ;;
  "")            ;;
  *) echo "usage: $0 [--python-only|--r-only]" >&2; exit 2 ;;
esac

say() { printf '\n==> %s\n' "$1"; }

# --------------------------------------------------------------------------- Python
if [ "$DO_PY" = 1 ]; then
  say "Python"
  PY_BIN="$(tools/detect_python.sh || true)"
  if [ -z "$PY_BIN" ]; then
    MIN="$(sed -n 's/^.*python_requires: *>= *\([0-9][0-9.]*\).*/\1/p' requirements.txt | head -1)"
    cat >&2 <<EOF
  No interpreter on PATH meets the required Python >= ${MIN:-3.12}.
  Install one and re-run, e.g.:
      brew install python@${MIN:-3.12}          # macOS
      sudo apt-get install python${MIN:-3.12}-venv   # Debian/Ubuntu
  Then:  tools/setup_env.sh
EOF
    exit 1
  fi
  echo "  interpreter: $PY_BIN ($("$PY_BIN" -V 2>&1))"

  # Recreated, not reused: a stale .venv can keep superseded packages, which the version guard catches.
  rm -rf "$VENV"
  "$PY_BIN" -m venv "$VENV"
  "$VENV/bin/python" -m pip install --quiet --upgrade pip
  "$VENV/bin/python" -m pip install --quiet -r requirements.txt
  echo "  installed into .venv/"
  "$VENV/bin/python" src/analysis/_version_guard.py
fi

# --------------------------------------------------------------------------- R
if [ "$DO_R" = 1 ]; then
  say "R prerequisites"

  if [ "$(uname -s)" = "Darwin" ]; then
  # Compiles trivial C++ to catch missing or shadowed standard headers, not merely a missing compiler.
  tmp="$(mktemp -t ignite_cxx).cpp"; printf '#include <cstdlib>\nint main(){return 0;}\n' > "$tmp"
  MAKEVARS=""
  if ! clang++ -c "$tmp" -o "${tmp%.cpp}.o" >/dev/null 2>&1; then
    # Locates and test-compiles SDK libc++ first, so the fallback is never offered where it would not work.
    SDK="$(xcrun --show-sdk-path 2>/dev/null || true)"
    FALLBACK_OK=0
    if [ -n "$SDK" ] && [ -f "$SDK/usr/include/c++/v1/cstdlib" ]; then
      clang++ -nostdinc++ -isystem "$SDK/usr/include/c++/v1" -isysroot "$SDK" \
              -c "$tmp" -o "${tmp%.cpp}.o" >/dev/null 2>&1 && FALLBACK_OK=1
    fi
    rm -f "$tmp" "${tmp%.cpp}.o"

    if [ "$FALLBACK_OK" = 1 ] && [ -n "${IGNITE_CXX_SDK_FALLBACK:-}" ] \
       && [ "${IGNITE_CXX_SDK_FALLBACK}" != "0" ]; then
      CXX_FALLBACK_INC="-nostdinc++ -isystem $SDK/usr/include/c++/v1 -isysroot $SDK"
      # OBJCXXFLAGS, not CXXFLAGS, must be set: systemfonts' font backend is Objective-C++.
      cat >&2 <<EOF
  ⚠ C++ toolchain is broken; building against the SDK's libc++ instead.
    IGNITE_CXX_SDK_FALLBACK is set, so setup is continuing with:
        $INC
    This is a WORKAROUND, not a fix. It changes how every R package here compiles, and your
    broken toolchain will still break other projects. Repair it when you can:
        sudo rm -rf /Library/Developer/CommandLineTools && sudo xcode-select --install
EOF
    else
      cat >&2 <<'EOF'
  C++ TOOLCHAIN BROKEN -- clang++ cannot compile #include <cstdlib>.

  Several pinned R packages are C++ and cannot build until this is fixed. The usual cause on
  macOS is a Command Line Tools install older than the OS, whose own c++/v1 directory is
  incomplete and shadows the SDK's complete copy.

  This needs your password, so it cannot be done for you:
      sudo rm -rf /Library/Developer/CommandLineTools && sudo xcode-select --install

  Then re-run: tools/setup_env.sh --r-only
EOF
      if [ "$FALLBACK_OK" = 1 ]; then
        cat >&2 <<'EOF'

  To proceed WITHOUT repairing it -- a workaround, not a fix, and it changes how every R
  package here compiles:
      IGNITE_CXX_SDK_FALLBACK=1 tools/setup_env.sh --r-only
EOF
      fi
      exit 1
    fi
  else
    rm -f "$tmp" "${tmp%.cpp}.o"
    echo "  C++ toolchain: OK"
  fi

  # See REPRODUCIBILITY.md, R (stage 2, `make figures`), for the freetype/harfbuzz probe details.
  BP="$(brew --prefix 2>/dev/null || echo /opt/homebrew)"
  probe_cflags="-I$BP/include -I$BP/include/freetype2 -I$BP/include/harfbuzz"
  if command -v pkg-config >/dev/null 2>&1; then
    probe_cflags="$probe_cflags $(pkg-config --cflags freetype2 harfbuzz 2>/dev/null || true)"
  fi
  missing=""
  for probe in "ft2build.h:freetype" "hb-ft.h:harfbuzz"; do
    hdr="${probe%%:*}"; brewname="${probe#*:}"
    printf '#include <%s>\nint main(){return 0;}\n' "$hdr" > "$tmp.c"
    cc -c "$tmp.c" -o "$tmp.o" $probe_cflags >/dev/null 2>&1 || missing="$missing $brewname"
  done
  rm -f "$tmp.c" "$tmp.o"
  if [ -n "$missing" ]; then
    cat >&2 <<EOF
  NATIVE HEADERS MISSING:$missing

  ragg and textshaping cannot build from source without them. Install (no password needed
  with Homebrew):
      brew install freetype libpng jpeg-turbo libtiff harfbuzz fribidi

  Then re-run: tools/setup_env.sh --r-only
EOF
    exit 1
  fi
  echo "  native headers (freetype, harfbuzz): OK"

  # Homebrew include/lib paths go to R regardless of the C++ fallback, since R's Makeconf ignores Homebrew and link fails otherwise.
  MAKEVARS="$(mktemp -t ignite_makevars)"
  {
    for std in "" 11 14 17 20; do
      echo "CXX${std}FLAGS = -falign-functions=64 -Wall -g -O2 ${CXX_FALLBACK_INC:-}"
    done
    echo "OBJCXXFLAGS = -g -O2 ${CXX_FALLBACK_INC:-}"
    echo "CPPFLAGS = -I$BP/include -I$BP/include/freetype2 -I$BP/include/harfbuzz"
    echo "LDFLAGS = -L$BP/lib"
  } >> "$MAKEVARS"

  else
  # Linux: system compiler and headers; package names are Debian/Ubuntu's.
  tmp="$(mktemp "${TMPDIR:-/tmp}/ignite_cxx.XXXXXX")"
  MAKEVARS=""
  printf '#include <cstdlib>\nint main(){return 0;}\n' > "$tmp.cpp"
  if ! "${CXX:-c++}" -c "$tmp.cpp" -o "$tmp.o" >/dev/null 2>&1; then
    rm -f "$tmp" "$tmp.cpp" "$tmp.o"
    cat >&2 <<'EOF'
  C++ TOOLCHAIN MISSING -- c++ cannot compile #include <cstdlib>.

  Several pinned R packages are C++. On Debian/Ubuntu:
      sudo apt-get install build-essential

  Then re-run: tools/setup_env.sh --r-only
EOF
    exit 1
  fi
  echo "  C++ toolchain: OK"

  probe_cflags=""
  if command -v pkg-config >/dev/null 2>&1; then
    probe_cflags="$(pkg-config --cflags freetype2 harfbuzz 2>/dev/null || true)"
  fi
  missing=""
  for hdr in ft2build.h hb-ft.h fribidi.h png.h jpeglib.h tiffio.h fontconfig/fontconfig.h; do
    # stdio.h first: jpeglib.h uses size_t and FILE without including them, so it cannot be probed alone.
    printf '#include <stdio.h>\n#include <%s>\nint main(){return 0;}\n' "$hdr" > "$tmp.c"
    cc -c "$tmp.c" -o "$tmp.o" $probe_cflags $(pkg-config --cflags fribidi 2>/dev/null || true) \
      >/dev/null 2>&1 || missing="$missing $hdr"
  done
  rm -f "$tmp" "$tmp.c" "$tmp.cpp" "$tmp.o"
  if [ -n "$missing" ]; then
    cat >&2 <<EOF
  NATIVE HEADERS MISSING:$missing

  ragg, textshaping and systemfonts cannot build from source without them. On Debian/Ubuntu:
      sudo apt-get install libfreetype6-dev libharfbuzz-dev libfribidi-dev libpng-dev \\
        libjpeg-dev libtiff-dev libfontconfig1-dev libcurl4-openssl-dev libxml2-dev

  Then re-run: tools/setup_env.sh --r-only
EOF
    exit 1
  fi
  echo "  native headers (freetype, harfbuzz, fribidi, png, jpeg, tiff, fontconfig): OK"

  fi

  say "R packages -> .rlib/"
  mkdir -p "$RLIB"
  # Installs binaries first, then source for shortfalls: CRAN's macOS binaries lag source releases, landing a patch behind the pins.
  R_LIBS_USER="$RLIB" R_MAKEVARS_USER="${MAKEVARS:-}" \
  PKG_CONFIG_PATH="$(brew --prefix 2>/dev/null || echo /opt/homebrew)/lib/pkgconfig:${PKG_CONFIG_PATH:-}" \
  IGNITE_RLIB="$RLIB" Rscript - <<'RS'
lib <- Sys.getenv("IGNITE_RLIB"); .libPaths(c(lib, .libPaths()))
repo <- "https://cloud.r-project.org"

# Parses the pins programmatically rather than restating them, using the same pkg==X.Y.Z format check_r_versions.R relies on.
raw   <- readLines("R-requirements.txt", warn = FALSE)
lines <- trimws(sub("#.*$", "", raw))
pins  <- lines[grepl("^[A-Za-z0-9._-]+==[0-9]", lines)]
want  <- setNames(sub("^[^=]+==", "", pins), sub("==.*$", "", pins))

short <- function() {
  bad <- character(0)
  for (p in names(want)) {
    have <- tryCatch(as.character(packageVersion(p, lib.loc = lib)), error = function(e) NA)
    if (is.na(have) || have != want[[p]]) bad <- c(bad, p)
  }
  bad
}

# Pass 1 installs binaries, which succeeds only when the pinned version matches the current CRAN release.
todo <- short()
if (length(todo)) install.packages(todo, lib = lib, repos = repo, quiet = TRUE)

# Pass 2 installs from source, since CRAN's macOS binaries can lag the source release by about a patch version.
todo <- short()
if (length(todo)) {
  cat("  binaries lag CRAN for:", paste(todo, collapse = ", "), "-- building from source\n")
  install.packages(todo, lib = lib, repos = repo, type = "source", quiet = TRUE)
}

# See REPRODUCIBILITY.md, R (stage 2, `make figures`), for why pass 3 is needed.
todo <- short()
if (length(todo)) {
  cat("  pinned below current CRAN for:", paste(todo, collapse = ", "), "-- fetching from the Archive\n")
  # Installs the tarball directly: remotes::install_version's lib argument includes the read-only R framework directory and fails.
  for (p in todo) {
    v <- want[[p]]
    for (u in c(sprintf("%s/src/contrib/Archive/%s/%s_%s.tar.gz", repo, p, p, v),
                sprintf("%s/src/contrib/%s_%s.tar.gz", repo, p, v))) {
      f <- file.path(tempdir(), basename(u))
      ok <- tryCatch({ download.file(u, f, quiet = TRUE); TRUE }, error = function(e) FALSE,
                     warning = function(w) FALSE)
      if (!ok) next
      try(install.packages(f, lib = lib, repos = NULL, type = "source", quiet = TRUE),
          silent = TRUE)
      break
    }
  }
}

todo <- short()
if (length(todo)) {
  cat("\n  COULD NOT REACH THE PINNED VERSION FOR:", paste(todo, collapse = ", "), "\n")
  for (p in todo) {
    have <- tryCatch(as.character(packageVersion(p, lib.loc = lib)), error = function(e) "ABSENT")
    cat(sprintf("    %-12s want %-8s got %s\n", p, want[[p]], have))
  }
  quit(status = 1)
}
cat("  all", length(want), "pinned R packages installed into .rlib/\n")
RS
fi

say "done"
cat <<'EOF'
  Both stacks are project-local: .venv/ and .rlib/, and both are gitignored.
  The Makefile picks them up automatically -- no PY= or R_LIBS_USER= needed.

  Check:  make check-versions
EOF
