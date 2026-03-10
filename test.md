```bash
# 1. Có file cài đặt do Sentieon cung cấp, ví dụ:
# sentieon-genomics-202503.02.tar.gz

mkdir -p $HOME/opt/sentieon
tar xzf sentieon-genomics-202503.02.tar.gz -C $HOME/opt/sentieon

# 2. Thêm binary vào PATH
export SENTIEON_BIN_DIR=$HOME/opt/sentieon/sentieon-genomics-202503.02/bin
export PATH="${SENTIEON_BIN_DIR}:$PATH"

# 3. Cấu hình license
export SENTIEON_LICENSE=/path/to/Hanoi_University_Of_Science_And_Technology_eval.lic
# hoặc:
# export SENTIEON_LICENSE=license-server-host:port

# 4. Kiểm tra
which sentieon
sentieon --help | head
which sentieon-cli
```