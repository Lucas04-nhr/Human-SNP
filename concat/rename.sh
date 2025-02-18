#!/bin/bash

# 遍历当前目录下的所有文件
for file in *; do
  # 检查是否是文件
  if [ -f "$file" ]; then
    # 使用正则表达式匹配 P 后面的数字
    if [[ $file =~ P([0-9]) ]]; then
      # 提取数字并补足为两位
      number=$(printf "%02d" "${BASH_REMATCH[1]}")
      new_file="${file/P${BASH_REMATCH[1]}/P$number}"
      
      # 如果新文件名与原文件名不同，执行重命名
      if [ "$file" != "$new_file" ]; then
        mv "$file" "$new_file"
        echo "Renamed '$file' to '$new_file'"
      fi
    fi
  fi
done
