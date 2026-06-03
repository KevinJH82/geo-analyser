#!/bin/bash
# 启动遥感数据预处理 Web 应用

cd "$(dirname "$0")"

echo "🚀 启动遥感数据预处理工作台..."
echo ""
echo "📍 Web 界面: http://127.0.0.1:5001"
echo "🛑 停止服务: Ctrl+C"
echo ""

python3 app.py
