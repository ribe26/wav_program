import re
import pandas as pd
from collections import Counter

# .txtファイルのパス
txt_file_path = 'datas.txt'

# テキストファイルを読み込む
with open(txt_file_path, 'r', encoding='utf-8') as file:
    data = file.read()

# データを分割
blocks = data.strip().split('-----------------------------')

# パラメータを格納するリストとCounterを準備
start_end_pairs = []
filter_lengths = []
min_fm_values = []
rt_values = []
range_labels = []  # minFs~maxFs の範囲ラベル

# 各ブロックからパラメータを抽出
for block in blocks:
    # startFs, endFs, FilterLength, minFm, RTを抽出
    startFs = re.search(r'startFs(\d+)', block)
    endFs = re.search(r'endFs(\d+)', block)
    filter_length = re.search(r'FilterLength(\d+)', block)
    minFm = re.search(r'minFm(\d+)', block)
    rt = re.search(r'RT(\d+(\.\d+)?)', block)  # 整数または少数の両方に対応する正規表現
    
    if startFs and endFs and filter_length and minFm and rt:
        # FilterLengthが1024の場合は除外
        filter_length_value = int(filter_length.group(1))
        if filter_length_value == 1024:
            continue
        
        # RTの値を整数または浮動小数点数として統一
        rt_value = float(rt.group(1))  # 整数も少数も float に統一
        
        # RTが1の場合は除外
        if rt_value == 1:
            continue
        
        # startFs と endFs のペアを追加
        start_end_pairs.append((int(startFs.group(1)), int(endFs.group(1))))
        filter_lengths.append(filter_length_value)
        min_fm_values.append(int(minFm.group(1)))
        rt_values.append(rt_value)
        
        # minFs~maxFs の範囲ラベルを作成
        range_labels.append(f"minFs{startFs.group(1)}~maxFs{endFs.group(1)}")

# 並べ替え: startFs と endFs のペアをソート
start_end_pairs.sort()

# 出現回数をカウント
filter_length_counts = Counter(filter_lengths)
min_fm_counts = Counter(min_fm_values)
rt_counts = Counter(rt_values)
range_counts = Counter(range_labels)  # minFs~maxFs の範囲のカウント

# DataFrameの作成
df = pd.DataFrame({
    'startFs': [pair[0] for pair in start_end_pairs],
    'endFs': [pair[1] for pair in start_end_pairs],
    'FilterLength': filter_lengths,
    'minFm': min_fm_values,
    'RT': rt_values,
    'Range': range_labels  # minFs~maxFsの範囲を追加
})

# 出現回数のデータフレームを作成
count_df = pd.DataFrame({
    'FilterLength Counts': filter_length_counts,
    'minFm Counts': min_fm_counts,
    'RT Counts': rt_counts
})

# minFs~maxFsの範囲のカウントを追加
range_count_df = pd.DataFrame({
    'Range Counts': range_counts
})

# エクセルファイルに書き込む
output_file_path = 'output_with_range_counts_excluding_RT1.xlsx'
with pd.ExcelWriter(output_file_path, engine='openpyxl') as writer:
    # 並べ替えたデータを1つのシートに保存
    df.to_excel(writer, sheet_name='Sorted Data', index=False)
    # 出現回数を別のシートに保存
    count_df.to_excel(writer, sheet_name='Counts', index=True)
    # minFs~maxFsの範囲のカウントを別のシートに保存
    range_count_df.to_excel(writer, sheet_name='Range Counts', index=True)

print(f"エクセルファイルに保存しました: {output_file_path}")
