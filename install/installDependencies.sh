# Install dependencies

# Check the operating system
if [[ "$(uname)" == "Darwin" ]]; then # macOS
    brew install wget git boost tbb cmake protobuf
else # Assuming Ubuntu otherwise
    sudo apt install -y wget git build-essential libboost-all-dev libtbb2 cmake protobuf-compiler
fi