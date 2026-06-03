"""
生成测试用的合成遥感数据 (.npy 格式)
用于快速演示和验证 Web 应用功能
"""

import numpy as np
from pathlib import Path

def generate_synthetic_image(
    bands: int = 7,
    rows: int = 256,
    cols: int = 256,
    scene_type: str = "mixed"
) -> np.ndarray:
    """
    生成合成多光谱影像

    Parameters
    ----------
    bands : int
        波段数
    rows, cols : int
        影像大小
    scene_type : str
        场景类型: "mixed" (混合), "forest" (森林), "urban" (城市)

    Returns
    -------
    np.ndarray
        DN 值影像 (bands, rows, cols), dtype=uint16
    """
    np.random.seed(42)
    image = np.random.uniform(100, 3000, (bands, rows, cols)).astype(np.uint16)

    if scene_type == "mixed":
        # 植被区块（高 NIR，低 Red）
        image[4, 50:100, 50:100] = 4500      # NIR
        image[3, 50:100, 50:100] = 1000      # Red

        # 水体区块（高 Green，低 NIR）
        image[2, 120:160, 120:160] = 3000    # Green
        image[4, 120:160, 120:160] = 500     # NIR

        # 城镇/建筑（高 SWIR，中 Red）
        image[3, 200:240, 200:240] = 3500    # Red
        image[4, 200:240, 200:240] = 2000    # NIR

    elif scene_type == "forest":
        # 主要是高 NDVI 区域
        image[4, 50:200, 50:200] = 5000      # 高 NIR
        image[3, 50:200, 50:200] = 800       # 低 Red
        image[2, 50:200, 50:200] = 1500      # 中等 Green

    elif scene_type == "urban":
        # 主要是建筑和道路
        image[5, :, :] = 2500                # 高 SWIR1
        image[4, :, :] = 1500                # 中等 NIR
        image[3, :, :] = 2800                # 高 Red

    return image


def create_demo_dataset():
    """创建演示数据集"""
    output_dir = Path("demo_data")
    output_dir.mkdir(exist_ok=True)

    scenes = [
        ("scene1", "mixed"),
        ("scene2", "forest"),
        ("scene3", "urban"),
    ]

    print("🎬 生成演示数据...")
    for scene_name, scene_type in scenes:
        image = generate_synthetic_image(scene_type=scene_type)
        file_path = output_dir / f"{scene_name}.npy"
        np.save(file_path, image)
        print(f"  ✓ {file_path} ({image.nbytes / 1024 / 1024:.1f}MB)")

    print(f"\n✅ 已生成 {len(scenes)} 个测试文件")
    print(f"📂 输入目录: {output_dir.absolute()}")
    print(f"\n📝 Web 应用中:")
    print(f"   输入目录输入: {output_dir.absolute()}")
    print(f"   输出目录输入: ./demo_output/")


if __name__ == "__main__":
    create_demo_dataset()
