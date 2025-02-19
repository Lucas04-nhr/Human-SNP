#!/bin/bash

# 遍历当前目录下的所有文件
for file in *; do
    # 检查文件是否存在
    if [ -f "$file" ]; then
        # 使用正则表达式匹配文件名中的模式
        if [[ "$file" =~ ^(.*\.P)([0-9]+)(\..*)$ ]]; then
            # 获取匹配的部分
            prefix="${BASH_REMATCH[1]}"
            number="${BASH_REMATCH[2]}"
            suffix="${BASH_REMATCH[3]}"
            
            # 检查数字是否为一位数
            if [ ${#number} -eq 1 ]; then
                # 将一位数格式化为两位数
                formatted_number=$(printf "%02d" "$number")
                
                # 构建新的文件名
                new_filename="${prefix}${formatted_number}${suffix}"
                
                # 执行重命名操作
                mv "$file" "$new_filename"
                
                # 打印重命名信息
                echo "Renamed: $file -> $new_filename"
            fi
        fi
    fi
done

echo "Batch renaming completed."

