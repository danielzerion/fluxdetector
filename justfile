# -*-Makefile-*-

set positional-arguments := true

test:
  just run --beam-on 10

build:
  meson setup build/app src
  meson compile -C build/app
  meson install -C build/app

run *ARGS: build
  #!/usr/bin/env sh
  sh execute-with-nixgl-if-needed.sh ./build/app/fluxdetector "$@"
  exit $?

clean:
  rm build -rf
