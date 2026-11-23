# Install dependencies

# Check the operating system
if [[ "$(uname)" == "Darwin" ]]; then # macOS
    brew install wget git cmake tbb
else # Assuming Ubuntu otherwise
    sudo apt update
    sudo apt install -y wget git build-essential cmake libtbb-dev zlib1g-dev
fi
