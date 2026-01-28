#!/bin/bash
# 批量将 SEED/stdlibs 下的结构转换为 VASP 计算目录（POSCAR/INCAR/POTCAR）
# 依赖：cabal（格式转换）、DFT/dyn_vasp_in、DFT/POTCAR-*，可选任务脚本（如 sbatch）

set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SRC_DIR="$SCRIPT_DIR/stdlibs"
DFT_DIR="$SCRIPT_DIR/../DFT"
CALC_TYPE="cell-relax"
TASK_SCRIPT=""
AUTO_YES=0
SUBMIT=0

show_help() {
    cat <<'EOF'
用法: prepare_stdlib_vasp.sh [选项]

批量读取 SEED/stdlibs 下的结构（*.res/*.vasp/*.poscar/*.cif），为每个文件创建同名工作目录，
生成 POSCAR/INCAR/POTCAR，并可附带任务脚本。

选项:
  -s, --source DIR        指定结构源目录（默认: SEED/stdlibs）
  -ct, --calculation TYPE 计算类型: scf 或 cell-relax（默认: cell-relax）
  -ts, --task SCRIPT      需拷贝的任务脚本（如 sbatch 提交脚本）
  --submit                准备完成后即提交任务（需配合 -ts 指定脚本）
  -y, --yes               跳过参数确认
  -h, --help              显示本帮助

行为说明:
  - 工作目录命名为“去后缀文件名”，已存在将直接删除后重建并给出 WARNING。
  - POSCAR 由 cabal 转换生成；INCAR 来源于 DFT/dyn_vasp_in；POTCAR 由 DFT/POTCAR-*
    按 POSCAR 第 6 行元素顺序拼接。缺少 POTCAR 时仅提示警告。
  - 支持文件类型：.res/.cif/.vasp/.poscar（大小写不敏感）。

示例:
  1) 默认源、cell-relax：  ./prepare_stdlib_vasp.sh
  2) 指定源目录与 SCF：    ./prepare_stdlib_vasp.sh -s /path/to/stdlibs -ct scf
  3) 拷贝任务脚本并跳过确认：
     ./prepare_stdlib_vasp.sh -ts ../DFT/sub.sh -y
EOF
}

err() { echo "[ERROR] $*" >&2; }
warn() { echo "[WARN] $*" >&2; }

convert_structure() {
    # $1: src, $2: dst (POSCAR)
    local src="$1" dst="$2" ext
    ext="$(echo "${src##*.}" | tr '[:upper:]' '[:lower:]')"
    case "$ext" in
        res)
            cabal res poscar <"$src" >"$dst"
            ;;
        cif)
            cabal cif poscar <"$src" >"$dst"
            ;;
        poscar|vasp)
            cp "$src" "$dst"
            ;;
        *)
            return 1
            ;;
    esac
}

build_potcar() {
    # $1: workdir, $2: poscar
    local workdir="$1" poscar="$2" types_line types upf=() t pot
    types_line=$(sed -n '6p' "$poscar" | tr -s '[:space:]' ' ')
    if [[ -z "$types_line" ]]; then
        warn "$workdir: POSCAR 第6行无法解析元素类型，跳过 POTCAR 生成"
        return 1
    fi
    read -r -a types <<<"$types_line"
    for t in "${types[@]}"; do
        pot="$DFT_DIR/POTCAR-$t"
        if [[ -f "$pot" ]]; then
            upf+=("$pot")
        else
            warn "$workdir: 未找到 $pot"
        fi
    done
    if [[ ${#upf[@]} -eq 0 ]]; then
        warn "$workdir: 无可用 POTCAR，跳过生成"
        return 1
    fi
    cat "${upf[@]}" >"$workdir/POTCAR"
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        -s|--source)
            SRC_DIR="$2"
            shift 2
            ;;
        -ct|--calculation)
            CALC_TYPE="$2"
            shift 2
            ;;
        -ts|--task)
            TASK_SCRIPT="$2"
            shift 2
            ;;
        --submit)
            SUBMIT=1
            shift
            ;;
        -y|--yes)
            AUTO_YES=1
            shift
            ;;
        -h|--help)
            show_help
            exit 0
            ;;
        *)
            err "未知参数: $1"
            show_help
            exit 1
            ;;
    esac
done

if [[ "$CALC_TYPE" != "scf" && "$CALC_TYPE" != "cell-relax" ]]; then
    err "计算类型必须为 scf 或 cell-relax"
    exit 1
fi

if [[ ! -x "$DFT_DIR/dyn_vasp_in" ]]; then
    err "未找到或不可执行: $DFT_DIR/dyn_vasp_in"
    exit 1
fi

if ! command -v cabal >/dev/null 2>&1; then
    err "缺少 cabal 命令，请先安装或加载环境"
    exit 1
fi

if [[ ! -d "$SRC_DIR" ]]; then
    err "源目录不存在: $SRC_DIR"
    exit 1
fi

mapfile -t INPUTS < <(find "$SRC_DIR" -maxdepth 1 -type f \( -iname "*.res" -o -iname "*.vasp" -o -iname "*.poscar" -o -iname "*.cif" \) | sort)

if [[ ${#INPUTS[@]} -eq 0 ]]; then
    err "在 $SRC_DIR 未找到可处理的结构文件"
    exit 1
fi

echo "SOURCE      : $SRC_DIR"
echo "WORK ROOT   : $SCRIPT_DIR"
echo "CALC TYPE   : $CALC_TYPE"
echo "DFT TEMPLATE: $DFT_DIR/dyn_vasp_in"
if [[ -n "$TASK_SCRIPT" ]]; then
    echo "TASK SCRIPT : $TASK_SCRIPT"
fi
echo "FILES       : ${#INPUTS[@]}"

if [[ $AUTO_YES -eq 0 ]]; then
    read -rp "确认继续并可能覆盖同名目录? [y/N] " u
    if [[ "$u" != "y" && "$u" != "Y" ]]; then
        echo "已取消"
        exit 0
    fi
fi

for src in "${INPUTS[@]}"; do
    bname="$(basename "$src")"
    work="${bname%.*}"
    workdir="$SCRIPT_DIR/$work"
    if [[ -d "$workdir" ]]; then
        warn "$workdir 已存在，执行覆盖"
    fi
    mkdir -p "$workdir"

    if ! convert_structure "$src" "$workdir/POSCAR"; then
        warn "$bname: 不支持的格式，跳过"
        continue
    fi

    "$DFT_DIR/dyn_vasp_in" -t "$CALC_TYPE" >"$workdir/INCAR"
    build_potcar "$workdir" "$workdir/POSCAR" || true

    if [[ -n "$TASK_SCRIPT" ]]; then
        cp "$TASK_SCRIPT" "$workdir"/
    fi

    if [[ $SUBMIT -eq 1 ]]; then
        if [[ -z "$TASK_SCRIPT" ]]; then
            warn "$workdir: 未指定任务脚本，无法提交"
        elif [[ ! -f "$workdir/$(basename "$TASK_SCRIPT")" ]]; then
            warn "$workdir: 找不到任务脚本 $(basename "$TASK_SCRIPT")，跳过提交"
        else
            (
                cd "$workdir"
                if command -v safe_sbatch >/dev/null 2>&1; then
                    safe_sbatch "$(basename "$TASK_SCRIPT")" || warn "$workdir: safe_sbatch 提交失败"
                else
                    warn "$workdir: 未找到 safe_sbatch，跳过提交"
                fi
            )
        fi
    fi

    echo "准备完成: $workdir"
done

echo "全部处理完成。"
