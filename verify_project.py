"""
整个项目的最终检查清单
"""

import os
import sys
from pathlib import Path

def check_project():
    """检查项目完整性"""

    base_dir = Path(".")

    print("="*60)
    print("🛰 遥感数据预处理工作台 - 项目检查")
    print("="*60)

    # 检查预处理模块
    print("\n📦 预处理模块:")
    modules = [
        ("atmospheric_correction.py", "大气校正"),
        ("geometric_correction.py", "几何校正"),
        ("interference_removal.py", "干扰剔除"),
    ]

    for module_file, desc in modules:
        path = base_dir / module_file
        if path.exists():
            size = path.stat().st_size / 1024
            print(f"  ✓ {module_file:<35} {desc:<15} ({size:.1f}KB)")
        else:
            print(f"  ✗ {module_file:<35} {desc:<15} (缺失)")

    # 检查 Web 应用
    print("\n🌐 Web 应用:")
    web_files = [
        ("app.py", "Flask 后端"),
        ("templates/index.html", "前端界面"),
        ("static/style.css", "样式表"),
    ]

    for file, desc in web_files:
        path = base_dir / file
        if path.exists():
            size = path.stat().st_size / 1024
            print(f"  ✓ {file:<35} {desc:<15} ({size:.1f}KB)")
        else:
            print(f"  ✗ {file:<35} {desc:<15} (缺失)")

    # 检查配置和文档
    print("\n📚 文档和配置:")
    doc_files = [
        ("README.md", "完整使用手册"),
        ("QUICKSTART.md", "快速开始指南"),
        ("requirements.txt", "Python 依赖"),
        ("run.sh", "启动脚本"),
        ("generate_demo_data.py", "演示数据生成器"),
    ]

    for file, desc in doc_files:
        path = base_dir / file
        if path.exists():
            size = path.stat().st_size / 1024
            print(f"  ✓ {file:<35} {desc:<15} ({size:.1f}KB)")
        else:
            print(f"  ✗ {file:<35} {desc:<15} (缺失)")

    # 检查演示数据
    print("\n🎬 演示数据:")
    demo_dir = base_dir / "demo_data"
    if demo_dir.exists():
        demo_files = list(demo_dir.glob("*.npy"))
        print(f"  ✓ demo_data/ 目录存在，包含 {len(demo_files)} 个文件")
        for f in demo_files:
            size = f.stat().st_size / 1024 / 1024
            print(f"    - {f.name} ({size:.1f}MB)")
    else:
        print(f"  ✗ demo_data/ 目录不存在（可运行 generate_demo_data.py 创建）")

    # 检查依赖
    print("\n📦 Python 依赖:")
    required = [
        ("flask", "Flask"),
        ("numpy", "NumPy"),
        ("scipy", "SciPy"),
        ("matplotlib", "Matplotlib"),
    ]

    for module, name in required:
        try:
            __import__(module)
            print(f"  ✓ {name:<20} 已安装")
        except ImportError:
            print(f"  ✗ {name:<20} 未安装 (运行: pip install {module})")

    # 总结
    print("\n" + "="*60)
    print("✅ 项目就绪！")
    print("="*60)

    print("\n📖 后续步骤:")
    print("  1. 运行: python3 generate_demo_data.py")
    print("  2. 运行: python3 app.py")
    print("  3. 打开: http://127.0.0.1:5000")
    print("  4. 查看 QUICKSTART.md 了解详细使用方法")
    print()

if __name__ == "__main__":
    check_project()
