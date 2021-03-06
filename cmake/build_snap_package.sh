#!/bin/bash
# OPTIONS:
#   $1 install directory

if [[ $# -ne 1 ]]; then
  echo "This script requires one argument: "
  echo "OPTIONS: "
  echo "  $1 install directory"
  exit
fi

sudo make install
mkdir -p build-snap && cd build-snap
mv ../snapcraft.yaml .
# create run command
mkdir -p bin
echo '#!/bin/bash' > bin/run
echo '# execute command given as first parameter' >> bin/run
echo '"$@"' >> bin/run
chmod +x bin/run

sudo docker pull snapcore/snapcraft
echo "Start building snap package"
sudo docker run --net=host -v $PWD:$PWD -v $1:$1 -e SNAPCRAFT_SETUP_CORE=1 -w $PWD snapcore/snapcraft snapcraft
