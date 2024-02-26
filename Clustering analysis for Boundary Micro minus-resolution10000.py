#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pybedtools
import pyBigWig
import os
import numpy as np
import pandas as pd


# In[2]:


#Create a function to filter valid HMR regions
def filter_valid_intervals(intervals, bw):
    valid_intervals = []
    chroms = bw.chroms()
    for interval in intervals:
        chrom, start, end = interval[0], int(interval[1]), int(interval[2])

        # Check if the chromosome exists in the bigwig file
        if chrom not in chroms:
            continue

        chrom_length = chroms[chrom]
        start = max(0, start)
        end = min(chrom_length, end)
        if start < end:
            valid_intervals.append((chrom, start, end))
    return valid_intervals


# In[3]:


HMR_file = "banana/data_hic/Input1_Boundary_microC/Boundary_minus_1000_resolution10000.bed"
TF_folder = "banana/data_hic/Input3_new/"
bigwig_file1 = "banana/data_hic/Input2_CTCT_Boundary/ChIPseq_GSE124667_MCF7_CTCF_NT_NT_rep0.single-bowtie2-hg38-raw-GR.100.bw"
bigwig_file2 = "banana/data_hic/Input2-Cohesin-ForBoundary/Rad21_ENCSR000BTQ_rep1.mpbl.100.bw"
h1d_file = "banana/data_hic/H1D_Micro_10000/minus_10000_CI.bedGraph"


# In[4]:


#convert TF files to pybedtools.BedTool format,
tf_files = [os.path.join(TF_folder, f) for f in os.listdir(TF_folder) if f.endswith('.bed')]
bed_files = [pybedtools.BedTool(f) for f in tf_files]


# In[5]:


# Filter HMR intervals
HMR = pybedtools.BedTool(HMR_file)
bw = pyBigWig.open(bigwig_file1)
filtered_HMR = filter_valid_intervals(HMR, bw)
filtered_HMR


# In[6]:


# Convert filtered_HMR list to BedTool object
filtered_HMR_bedtool = pybedtools.BedTool(filtered_HMR)


# In[7]:


# Calculate overlap_matrix
overlap_matrix = []
for tf in bed_files:
    overlaps = filtered_HMR_bedtool.intersect(tf, c=True)
    overlap_counts = [int(interval[-1]) for interval in overlaps]
    overlap_matrix.append(overlap_counts)


# In[8]:


overlap_matrix = np.array(overlap_matrix).T
overlap_matrix = (overlap_matrix - overlap_matrix.min()) / (overlap_matrix.max() - overlap_matrix.min())


# In[9]:


# Calculate response_vector
def calculate_avg_signal(bigwig_file, filtered_HMR):
    bw = pyBigWig.open(bigwig_file)
    response_vector = []

    for interval in filtered_HMR:
        chrom, start, end = interval

        # Make sure start and end positions are within chromosome bounds
        chrom_length = bw.chroms()[chrom]
        start = max(0, start)
        end = min(chrom_length, end)

        if start < end:
            avg_signal = bw.stats(chrom, start, end, type='mean')[0]
            if avg_signal is None:
                avg_signal = 0
        else:
            avg_signal = 0

        response_vector.append(avg_signal)

    bw.close()

    return response_vector

# Calculate response_vector
response_vector1 = calculate_avg_signal(bigwig_file1, filtered_HMR)
response_vector2 = calculate_avg_signal(bigwig_file2, filtered_HMR)


# In[10]:


#log transformation
response_vector1 = np.array(response_vector1) + 1
response_vector1 = np.log(response_vector1)
response_vector1 = (response_vector1 - response_vector1.min()) / (response_vector1.max() - response_vector1.min())

response_vector2 = np.array(response_vector2) + 1
response_vector2 = np.log(response_vector2)
response_vector2 = (response_vector2 - response_vector2.min()) / (response_vector2.max() - response_vector2.min())


# In[11]:


def calculate_avg_signal_h1d(h1d_file, filtered_HMR):
    bed_intervals = [line.strip().split() for line in open(h1d_file)]
    bed_dict = {}
    for bed_interval in bed_intervals:
        chrom, start, end = bed_interval[0], int(bed_interval[1]), int(bed_interval[2])
        value = float(bed_interval[3]) if len(bed_interval) == 4 else 0.0
        if chrom not in bed_dict:
            bed_dict[chrom] = []
        bed_dict[chrom].append((start, end, value))

    response_vector = []
    for interval in filtered_HMR:
        chrom, start, end = interval
        if chrom not in bed_dict:
            response_vector.append(0)
            continue

        signals = []
        weights = []
        for bed_start, bed_end, value in bed_dict[chrom]:
            if bed_start < end and bed_end > start:  # If overlapping
                overlap_start = max(start, bed_start)
                overlap_end = min(end, bed_end)
                overlap_length = overlap_end - overlap_start
                weights.append(overlap_length)
                signals.append(value)

        if len(signals) == 1:
            response_vector.append(signals[0])
        elif len(signals) == 2:
            weighted_avg_signal = (signals[0]*weights[0] + signals[1]*weights[1]) / (weights[0] + weights[1])
            response_vector.append(weighted_avg_signal)
        else:
            response_vector.append(0)

    return response_vector

# Calculate response_vector
response_vector_h1d = calculate_avg_signal_h1d(h1d_file, filtered_HMR)

# 确保 response_vector_h1d 是一个 NumPy 数组
response_vector_h1d = np.array(response_vector_h1d)


# 现在可以计算最小值和最大值，并进行归一化

response_vector_h1d = np.array(response_vector_h1d) + 1
response_vector_h1d = np.log(response_vector_h1d)
response_vector_h1d = (response_vector_h1d - response_vector_h1d.min()) / (response_vector_h1d.max() - response_vector_h1d.min())


# In[12]:


tf_files_filtered = [name for name in tf_files if isinstance(name, str)]
tf_files_simplified = [name.replace('banana/data_hic/Input3_new/', '').replace('.bed', '') for name in tf_files_filtered]

tf_files_simplified


# In[13]:


import numpy as np
from scipy.spatial.distance import cosine, euclidean
from scipy.stats import pearsonr
import pandas as pd

# 确保 overlap_matrix 是 NumPy 数组
# 如果它原本不是，可以使用下面的代码行转换它
# overlap_matrix = np.array(overlap_matrix)

num_columns = overlap_matrix.shape[1]

# 初始化列表来存储结果
results = []

# 计算度量并存储结果
for i in range(num_columns):
    for j in range(i+1, num_columns):
        col_i = overlap_matrix[:, i]
        col_j = overlap_matrix[:, j]

        # 计算各种度量
        cosine_sim = 1 - cosine(col_i, col_j)  # 计算余弦相似度

        # 添加到结果列表
        results.append({
            'Factor_1': tf_files_simplified[i],
            'Factor_2': tf_files_simplified[j],
            'Cosine_Similarity': cosine_sim,
        })

# 创建DataFrame
df_A = pd.DataFrame(results)

# 将DataFrame保存为CSV文件
df_A.to_csv('comparison_results.csv', index=False)

# 显示DataFrame
print(df_A)


# In[14]:


df_A = df_A.dropna()
df_A


# In[15]:


import numpy as np
from scipy.spatial.distance import cosine, cityblock
import pandas as pd

# 假设 overlap_matrix 是您的矩阵，tf_files_simplified 是因子名称列表
# 假设有多个响应向量
# response_vector1 = [...], response_vector2 = [...], response_vector_h1d = [...]

response_vectors = {
    'signal_vector1': response_vector1,
    'signal_vector2': response_vector2,
    'h1d_vector': response_vector_h1d
}

num_columns = overlap_matrix.shape[1]
factor_names = tf_files_simplified  # 确保这个列表的长度与 overlap_matrix 的列数相同

# 初始化 DataFrame
df_B = pd.DataFrame({'Factor': factor_names})

# 计算每个响应向量与每一列的余弦相似度
for column_name, response_vector in response_vectors.items():
    cosine_similarities = []
    manhattan_distances = []
    for i in range(overlap_matrix.shape[1]):
        factor_vector = overlap_matrix[:, i]
        cosine_sim = 1 - cosine(factor_vector, response_vector)
        cosine_similarities.append(cosine_sim)
        
        # Calculate Manhattan distance
        manh_dist = cityblock(factor_vector, response_vector)
        manhattan_distances.append(manh_dist)
    
    # Normalize Manhattan distances
    max_manh_dist = max(manhattan_distances)
    min_manh_dist = min(manhattan_distances)
    normalized_manhattan_distances = [1 - (manh_dist - min_manh_dist) / (max_manh_dist - min_manh_dist) if max_manh_dist != min_manh_dist else 0 for manh_dist in manhattan_distances]
    
    df_B[column_name + '_cosine'] = cosine_similarities
    df_B[column_name + '_manhattan'] = normalized_manhattan_distances

# Display DataFrame
print(df_B)


# In[16]:


df_B


# In[17]:


df_B = df_B.dropna()
# If you want to display all rows of the dataframe, you can set the display.max_rows option in pandas to None
pd.set_option('display.max_rows', None)

# Now when you print df_B, it should display all the rows.
df_B


# In[18]:


import numpy as np

def generate_random_vectors(num_vectors, vector_length):
    """生成随机向量"""
    return np.random.rand(num_vectors, vector_length)

def cosine_similarity(v1, v2):
    """计算两个向量之间的余弦相似性"""
    norm_v1 = np.linalg.norm(v1)
    norm_v2 = np.linalg.norm(v2)
    return np.dot(v1, v2) / (norm_v1 * norm_v2) if norm_v1 != 0 and norm_v2 != 0 else 0

def calculate_max_avg_cosine_similarity(overlap_matrix, num_random_vectors, vector_length):
    """
    对于每个因子，计算它与一系列随机向量之间的余弦相似性的最大值，然后计算这些最大值的平均值。
    """
    random_vectors = generate_random_vectors(num_random_vectors, vector_length)
    max_similarities = []

    for factor in overlap_matrix.T:  # 遍历每个因子
        similarities = [cosine_similarity(factor, random_vector) for random_vector in random_vectors]
        max_similarities.append(max(similarities))

    return np.mean(max_similarities)

# 示例用法
num_random_vectors = 100  # 随机向量的数量
vector_length = overlap_matrix.shape[0]  # 向量长度等于overlap_matrix的行数

cutoff_value = calculate_max_avg_cosine_similarity(overlap_matrix, num_random_vectors, vector_length)
print("Cutoff Value:", cutoff_value)


# In[19]:


import pandas as pd
import networkx as nx
from stellargraph import StellarGraph
from stellargraph.layer import GCN
from stellargraph.mapper import FullBatchNodeGenerator
from tensorflow.keras.models import Model
from tensorflow.keras import optimizers
import numpy as np
from sklearn.cluster import KMeans

# 创建网络
G = nx.Graph()

# 因子自身特征的列名
factor_features = ['signal_vector1_cosine', 'signal_vector1_manhattan', 'signal_vector2_cosine', 'signal_vector2_manhattan', 'h1d_vector_cosine', 'h1d_vector_manhattan']

# 初始化DataFrame以存储节点特征
node_data = pd.DataFrame()

# 阈值
threshold = cutoff_value

# 添加节点和边
for index, row in df_A.iterrows():
    # 添加节点
    G.add_node(row['Factor_1'])
    G.add_node(row['Factor_2'])

    # 只有当权重大于或等于阈值时，才添加边
    if row['Cosine_Similarity'] >= threshold:
        G.add_edge(row['Factor_1'], row['Factor_2'], weight=row['Cosine_Similarity'])

for index, row in df_B.iterrows():
    # 聚合节点特征到DataFrame
    node_name = row['Factor']  # 或者 row['Factor_2']
    features = row[factor_features]  # 假设我们只取factor1_features
    node_data = node_data.append(pd.Series(features, name=node_name))
# 去除可能的重复行，并填充缺失值

#node_data = node_data[~node_data.index.duplicated(keep='first')].fillna(0)

# 创建StellarGraph对象，使用node_data作为节点特征
stellar_G = StellarGraph.from_networkx(G, node_features=node_data)





# In[20]:


import numpy as np
from tensorflow.keras import optimizers, Model
from tensorflow.keras.callbacks import EarlyStopping

# 设置 GCN 模型
gcn_layer_sizes = [16, 16]
gcn_generator = FullBatchNodeGenerator(stellar_G, method="gcn")
gcn_model = GCN(
    layer_sizes=gcn_layer_sizes, 
    generator=gcn_generator, 
    activations=["relu", "relu"],
    dropout=0.0
)

# 创建 Keras 模型
x_inp, x_out = gcn_model.in_out_tensors()
embedding_model = Model(inputs=x_inp, outputs=x_out)

# 获取节点 ID
node_ids = node_data.index

# 创建伪目标和数据生成器
dummy_targets = np.zeros((len(node_ids), gcn_layer_sizes[-1]))
node_gen = gcn_generator.flow(node_ids, dummy_targets)

# 编译模型
embedding_model.compile(optimizer=optimizers.Adam(learning_rate=0.001), loss='mse')

# 设置早停法参数
early_stopping_monitor = EarlyStopping(
    monitor='loss',  # 监控训练损失
    patience=1,  # 在 3 个 epochs 内如果损失没有减少则停止训练
    restore_best_weights=True  # 恢复到最佳模型权重
)

# 训练模型，添加早停法回调
embedding_model.fit(
    node_gen,
    epochs=20,
    callbacks=[early_stopping_monitor]
)

# 提取节点嵌入
node_embeddings = embedding_model.predict(node_gen)


# In[21]:


node_embeddings


# In[22]:


from scipy.cluster.hierarchy import dendrogram, linkage
from matplotlib import pyplot as plt

# 提取节点嵌入
node_embeddings = embedding_model.predict(node_gen)

# 重塑嵌入数据以确保其为二维
node_embeddings_reshaped = node_embeddings.squeeze()

# 使用分层聚类进行链接
#linked = linkage(node_embeddings_reshaped, 'ward')
#np.save('linked_matrix_boundary_minus_resolution10000.npy', linked)
linked = np.load('linked_matrix_boundary_minus_resolution10000.npy')

# 准备因子名称列表
# 确保这个列表的顺序与 node_embeddings_reshaped 中的行顺序相同
factor_names = node_data.index

# 绘制树状图
plt.figure(figsize=(50, 30))
dendrogram(linked, 
           orientation='top', 
           distance_sort='descending', 
           labels=factor_names,  # 使用因子名称作为标签
           show_leaf_counts=True)
plt.title('Hierarchical Clustering Dendrogram')
plt.xlabel('Factor Name')
plt.ylabel('Distance')
plt.xticks(rotation=90, fontsize=20)  # 旋转标签以便于阅读

plt.show()


# In[23]:


import numpy as np
from scipy.cluster.hierarchy import dendrogram
from matplotlib import pyplot as plt

# Load your linkage matrix
linked = np.load('linked_matrix_boundary_minus_resolution10000.npy')

# Replace with your actual factor names
factor_names = node_data.index

# Plotting the dendrogram rotated by 90 degrees to improve readability of factor names
plt.figure(figsize=(65, 108))  # Adjust the size as needed
dendrogram(linked,
           orientation='left',  # Rotate the dendrogram to display horizontally
           labels=factor_names,  # Use factor names as labels
           show_leaf_counts=True)
plt.title('Hierarchical Clustering Dendrogram (rotated)')
plt.xlabel('Distance')
plt.ylabel('Factor Name')
plt.xticks(fontsize=10)  # Adjust the font size as needed for distance ticks
plt.yticks(fontsize=50)  # Adjust the font size as needed for factor names
plt.show()


# In[24]:


from sklearn.metrics import silhouette_score

# 使用 fcluster 从 linkage matrix 中提取聚类标签
from scipy.cluster.hierarchy import fcluster
cluster_labels = fcluster(linked, t=0.05, criterion='distance')

# 计算轮廓系数
silhouette_avg = silhouette_score(node_embeddings_reshaped, cluster_labels)
print("The average silhouette_score is :", silhouette_avg)


# In[25]:


from scipy.cluster.hierarchy import fcluster
import pandas as pd

# 假设你已经有了linkage matrix 'linked' 和因子名称列表 'factor_names'
# 设置一个阈值来定义簇的距离，这个值应该基于dendrogram
# 例如，使用树状图中较大距离的值来确定切割点，这里假设为25
cut_off_distance = 0.4

# 使用fcluster根据选定的距离进行分簇
cluster_labels = fcluster(linked, cut_off_distance, criterion='distance')

# 将簇标签与因子名称映射成一个字典
clusters = {}
for label, factor in zip(cluster_labels, factor_names):
    if label not in clusters:
        clusters[label] = []
    clusters[label].append(factor)

# 创建一个DataFrame来展示每个簇的因子列表
cluster_df = pd.DataFrame.from_dict(clusters, orient='index').transpose()

# 输出DataFrame的前几行
print(cluster_df)

# 如果你想将DataFrame保存为CSV文件
cluster_df.to_csv("cluster_factors_boundary_minus_resolution10000.csv", index=False)


# In[26]:


cluster_3_factors = cluster_df[1].dropna().tolist()

cluster_3_factors 


# In[27]:


# Prepare a DataFrame to store the average values for each cluster
cluster_averages = pd.DataFrame(columns=['Cluster', 'Avg_cosine_similarity_vector1', 'Avg_cosine_similarity_vector2', 'Avg_cosine_similarity_h1d'])

for cluster_id, factors in cluster_df.iteritems():
    # Select the rows corresponding to the current cluster's factors
    cluster_data = df_B[df_B['Factor'].isin(factors.dropna().tolist())]
    
    # Calculate the mean values for the cosine similarity vectors
    mean_values = cluster_data[['signal_vector1_cosine', 'signal_vector1_manhattan', 'signal_vector2_cosine', 'signal_vector2_manhattan', 'h1d_vector_cosine', 'h1d_vector_manhattan']].mean()
    
    # Append the results to the cluster_averages DataFrame
    cluster_averages = cluster_averages.append({
        'Cluster': cluster_id,
        'Avg_cosine_similarity_vector1': mean_values['signal_vector1_cosine'],
        'Avg_manhattan_vector1': mean_values['signal_vector1_manhattan'],
        'Avg_cosine_similarity_vector2': mean_values['signal_vector2_cosine'],
        'Avg_manhattan_vector2': mean_values['signal_vector2_manhattan'],
        'Avg_cosine_similarity_h1d': mean_values['h1d_vector_cosine'],
        'Avg_manhattan_h1d': mean_values['h1d_vector_manhattan']
    }, ignore_index=True)

# Output the results
print(cluster_averages)

cluster_averages.to_csv("cluster_averages_boundary_minus_resolution10000.csv", index=False)


# In[28]:


import pandas as pd

# 假设 df_B 已经存在，并且包含所有因子的数据
# 以下代码应替换为您实际的数据加载/定义代码
# df_B = pd.read_csv('path_to_your_data.csv')  # 如果您的数据是从CSV文件加载的话

# 选择包含cosine similarity 和 manhattan distance 的列
columns_of_interest = ['signal_vector1_cosine', 'signal_vector1_manhattan',
                       'signal_vector2_cosine', 'signal_vector2_manhattan', 
                       'h1d_vector_cosine', 'h1d_vector_manhattan']

# 计算每个因子的平均值
df_B['Average'] = df_B[columns_of_interest].mean(axis=1)

# 从小到大排序因子的平均值
sorted_factors = df_B.sort_values(by='Average')

# 输出结果
print(sorted_factors[['Factor', 'Average']])

# 保存到CSV文件
sorted_factors.to_csv("sorted_factors_averages_resolution10000.csv", index=False)


# In[29]:


# 过滤 df_B 以只包含 cluster_3_factors 中的因子
# 假设 cluster_3_factors 是一个包含因子名称的列表
filtered_df_B = df_B[df_B['Factor'].isin(cluster_3_factors)]

# 计算这些因子的平均值
filtered_df_B['Average'] = filtered_df_B[columns_of_interest].mean(axis=1)

# 从小到大排序这些因子的平均值
sorted_filtered_factors = filtered_df_B.sort_values(by='Average')

# 输出结果
print(sorted_filtered_factors[['Factor', 'Average']])


# In[30]:


import seaborn as sns
import matplotlib.pyplot as plt

# 绘制箱型图
plt.figure(figsize=(8, 6)) # 根据需要调整图形大小
sns.boxplot(y=sorted_filtered_factors['Average'])

# 可选：添加更多定制化设置
plt.title('Boxplot of Average Values in Cluster2') # 添加标题
plt.xlabel('Cluster2') # 设置x轴标签
plt.ylabel('Average Value') # 设置y轴标签

# 显示图形
plt.show()


# In[31]:


# 修改后的代码
# Prepare a DataFrame to store the average values for each cluster
cluster_averages = pd.DataFrame(columns=['Cluster', 
                                         'Avg_cosine_similarity_vector1', 
                                         'Avg_manhattan_vector1',
                                         'Avg_cosine_similarity_vector2',
                                         'Avg_manhattan_vector2',
                                         'Avg_cosine_similarity_h1d',
                                         'Avg_manhattan_h1d',
                                         'Overall_Avg'])

for cluster_id, factors in cluster_df.iteritems():
    # Select the rows corresponding to the current cluster's factors
    cluster_data = df_B[df_B['Factor'].isin(factors.dropna().tolist())]
    
    # Calculate the mean values for the cosine similarity vectors
    mean_values = cluster_data[['signal_vector1_cosine', 
                                'signal_vector1_manhattan', 
                                'signal_vector2_cosine', 
                                'signal_vector2_manhattan', 
                                'h1d_vector_cosine', 
                                'h1d_vector_manhattan']].mean()
    
    # Calculate the overall average of the features for the cluster
    overall_avg = mean_values.mean()
    
    # Append the results to the cluster_averages DataFrame
    cluster_averages = cluster_averages.append({
        'Cluster': cluster_id,
        'Avg_cosine_similarity_vector1': mean_values['signal_vector1_cosine'],
        'Avg_manhattan_vector1': mean_values['signal_vector1_manhattan'],
        'Avg_cosine_similarity_vector2': mean_values['signal_vector2_cosine'],
        'Avg_manhattan_vector2': mean_values['signal_vector2_manhattan'],
        'Avg_cosine_similarity_h1d': mean_values['h1d_vector_cosine'],
        'Avg_manhattan_h1d': mean_values['h1d_vector_manhattan'],
        'Overall_Avg': overall_avg
    }, ignore_index=True)

# Output the results
print(cluster_averages)

# Save the DataFrame to a CSV file
cluster_averages.to_csv("cluster_averages_boundary_minus_resolution10000.csv", index=False)


# In[32]:


# 导入所需库
import numpy as np
import pandas as pd

# 假设 overlap_matrix 和 tf_files_simplified 已经被定义
# 假设 response_vector1, response_vector2, response_vector_h1d 已经被定义
# 假设 cluster_3_factors 已经包含了我们感兴趣的因子
# 假设 filtered_HMR 已经被定义

# 找出cluster_3_factors对应的overlap_matrix中的列
cluster_3_indices = [tf_files_simplified.index(factor) for factor in cluster_3_factors if factor in tf_files_simplified]
cluster_3_matrix = overlap_matrix[:, cluster_3_indices]

# 计算分子：overlap_matrix中相关列的平均值
numerator = np.sum(cluster_3_matrix, axis=1)

# 计算分母：三个response vectors的平均值
denominator = np.mean([response_vector1, response_vector2], axis=0)

# 计算比例
#ratios = (numerator + 1e-9) / (denominator + 1e-9)
ratios = (numerator + 1e-9) / 1

# 降序排序比例值并获取排序后的索引
sorted_indices_desc = np.argsort(ratios)[::-1]

# 使用排序后的索引找出filtered_HMR中对应的信息
sorted_filtered_HMR_desc = [filtered_HMR[i] for i in sorted_indices_desc]

# 创建一个DataFrame来更好地展示这些信息
sorted_filtered_HMR_df_desc = pd.DataFrame(sorted_filtered_HMR_desc, columns=['Chromosome', 'Start', 'End'])

# 将比例值也加入到DataFrame中，并按照降序排列
sorted_filtered_HMR_df_desc['Ratio'] = ratios[sorted_indices_desc]

# 显示或导出DataFrame
pd.set_option('display.max_rows', None)
print(sorted_filtered_HMR_df_desc)

# 可选：保存为CSV文件
#sorted_filtered_HMR_df_desc.to_csv("sorted_filtered_HMR_with_ratios_desc.csv", index=False)

# 显示 DataFrame 的前 100 行
print(sorted_filtered_HMR_df_desc.head(100))


# In[33]:


cluster_3_matrix.shape


# In[34]:


numerator.shape


# In[35]:


denominator.shape


# In[36]:


print(sorted_filtered_HMR_df_desc.head(100))


# In[37]:


import matplotlib.pyplot as plt

def plot_histogram(data, title, bins=50):
    """
    绘制数据的直方图。

    :param data: 数据数组。
    :param title: 直方图的标题。
    :param bins: 直方图的区间数量。
    """
    plt.figure(figsize=(10, 6))
    plt.hist(data, bins=bins, color='blue', alpha=0.7)
    plt.title(title)
    plt.xlabel('Value')
    plt.ylabel('Frequency')
    plt.grid(True)
    plt.show()

# 找出cluster_3_factors对应的overlap_matrix中的列
cluster_3_indices = [tf_files_simplified.index(factor) for factor in cluster_3_factors if factor in tf_files_simplified]
cluster_3_matrix = overlap_matrix[:, cluster_3_indices]

# 找到非零行的索引
non_zero_indices = np.where(cluster_3_matrix.any(axis=1))[0]

# 从response_vector1和response_vector2中提取对应的值
selected_values_vector1 = response_vector1[non_zero_indices]
selected_values_vector2 = response_vector2[non_zero_indices]

# 绘制直方图
plot_histogram(selected_values_vector1, 'Histogram of Selected BigWig Signal Vector 1')
plot_histogram(selected_values_vector2, 'Histogram of Selected BigWig Signal Vector 2')


# In[38]:


import matplotlib.pyplot as plt
import numpy as np

def plot_histogram(data, title, bins=50):
    """
    绘制数据的直方图。

    :param data: 数据数组。
    :param title: 直方图的标题。
    :param bins: 直方图的区间数量。
    """
    plt.figure(figsize=(10, 6))
    plt.hist(data, bins=bins, color='blue', alpha=0.7)
    plt.title(title)
    plt.xlabel('Log-transformed value')
    plt.ylabel('Frequency')
    plt.grid(True)
    plt.show()

# 为cluster_3_factors中的每个因子绘制直方图
for factor in cluster_3_factors:
    if factor in tf_files_simplified:
        # 找出因子对应的overlap_matrix中的列
        factor_index = tf_files_simplified.index(factor)
        factor_column = overlap_matrix[:, factor_index]

        # 找到非零行的索引
        non_zero_indices = np.where(factor_column > 0)[0]

        # 从response_vector1和response_vector2中提取对应的值，并计算平均值
        selected_values_vector1 = response_vector1[non_zero_indices]
        selected_values_vector2 = response_vector2[non_zero_indices]
        average_values = (selected_values_vector1 + selected_values_vector2) / 2


        # 绘制经过对数变换的值的直方图
        plot_histogram(average_values, f'Histogram of Log-transformed Average BigWig Signal for {factor}')



# In[39]:


# 假设 average_values 是我们之前计算的每个因子的平均信号值
for factor in cluster_3_factors:
    if factor in tf_files_simplified:
        factor_index = tf_files_simplified.index(factor)
        factor_column = overlap_matrix[:, factor_index]
        non_zero_indices = np.where(factor_column > 0)[0]
        
        selected_values = (response_vector1[non_zero_indices] + response_vector2[non_zero_indices]) / 2
        
        # 
        if len(selected_values) > 0:  # 确保数组不为空
            threshold_Negative = 0.4
            
            
            # 计算小于和大于阈值的区域数量
            low_count = np.sum(selected_values < threshold_Negative)
            
            print(f"Factor: {factor}")
            print(f"Regions below threshold: {low_count}")
            #print(f"Regions above threshold: {high_count}")
            print()


# In[40]:


import numpy as np

results = {}  # 用于存储结果的字典

for factor in cluster_3_factors:
    if factor in tf_files_simplified:
        factor_index = tf_files_simplified.index(factor)
        factor_column = overlap_matrix[:, factor_index]
        non_zero_indices = np.where(factor_column > 0)[0]
        selected_values = (response_vector1[non_zero_indices] + response_vector2[non_zero_indices]) / 2
        
        if len(selected_values) > 0:  # 确保数组不为空
            # 使用百分位数作为阈值
            threshold = np.percentile(selected_values, 30)  # 选取高于80%的值（相当于低于20%的阈值）

            # 根据阈值分类
            above_threshold = non_zero_indices[selected_values < threshold]

            # 获取区域的具体位置
            regions_above = [filtered_HMR[i] for i in above_threshold]

            # 存储结果
            results[factor] = {
                'threshold': threshold,
                'above_threshold': regions_above
            }

            print(f"Factor: {factor}")
            print(f"Regions above threshold: {len(regions_above)}")
            print()

# 如果你想查看具体结果，可以直接访问results字典
# 例如，查看某个因子的结果：
for factor, data in results.items():
    print(f"Results for {factor}:")
    print(f"Regions above threshold: {data['above_threshold']}")
    print()



# # Postive factor region detection

# In[41]:


cluster_4_factors = cluster_df[2].dropna().tolist()

cluster_4_factors 


# In[42]:


# 过滤 df_B 以只包含 cluster_3_factors 中的因子
# 假设 cluster_3_factors 是一个包含因子名称的列表
filtered_df_B = df_B[df_B['Factor'].isin(cluster_4_factors)]

# 计算这些因子的平均值
filtered_df_B['Average'] = filtered_df_B[columns_of_interest].mean(axis=1)

# 从小到大排序这些因子的平均值
sorted_filtered_factors1 = filtered_df_B.sort_values(by='Average')

# 输出结果
print(sorted_filtered_factors1[['Factor', 'Average']])


# In[54]:


import seaborn as sns
import matplotlib.pyplot as plt

# 假设 sorted_filtered_factors1 和 sorted_filtered_factors 是已经准备好的DataFrame

# 为了绘制箱型图，我们首先需要将两个DataFrame合并为一个，且包含一个类别标签列
sorted_filtered_factors1['Cluster'] = 'Cluster1'
sorted_filtered_factors['Cluster'] = 'Cluster2'
combined_data = pd.concat([sorted_filtered_factors1, sorted_filtered_factors], ignore_index=True)

# 绘制箱型图
plt.figure(figsize=(10, 8))
sns.boxplot(x='Cluster', y='Average', data=combined_data)

# 设置字体大小和标题
#plt.title("Feature’s average of transcription factors", fontsize=20)
plt.xlabel('', fontsize=30)  # 留空因为横坐标的 'Cluster1' 和 'Cluster2' 已经包含在类别里了
plt.ylabel('Feature’s average of transcription factors', fontsize=23)

# 调整横坐标的字体大小
plt.xticks(fontsize=30)

# 显示图形
plt.show()


# In[55]:


import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd  # 确保导入了 pandas

# 假设 sorted_filtered_factors1 和 sorted_filtered_factors 是已经准备好的DataFrame

# 为了绘制箱型图，我们首先需要将两个DataFrame合并为一个，且包含一个类别标签列
sorted_filtered_factors1['Cluster'] = 'Cluster1'
sorted_filtered_factors['Cluster'] = 'Cluster2'
combined_data = pd.concat([sorted_filtered_factors1, sorted_filtered_factors], ignore_index=True)

# 绘制箱型图
plt.figure(figsize=(10, 8))
sns.boxplot(x='Cluster', y='Average', data=combined_data,
            palette={'Cluster1': 'lightgreen', 'Cluster2': 'lightcoral'})  # 设置颜色

# 设置字体大小和标题
plt.xlabel('', fontsize=30)  # 留空因为横坐标的 'Cluster1' 和 'Cluster2' 已经包含在类别里了
plt.ylabel('Feature’s average of transcription factors', fontsize=23)

# 调整横坐标的字体大小
plt.xticks(fontsize=30)

# 显示图形
plt.show()


# In[44]:


import matplotlib.pyplot as plt

def plot_histogram(data, title, bins=50):
    """
    绘制数据的直方图。

    :param data: 数据数组。
    :param title: 直方图的标题。
    :param bins: 直方图的区间数量。
    """
    plt.figure(figsize=(10, 6))
    plt.hist(data, bins=bins, color='blue', alpha=0.7)
    plt.title(title)
    plt.xlabel('Value')
    plt.ylabel('Frequency')
    plt.grid(True)
    plt.show()

# 找出cluster_4_factors对应的overlap_matrix中的列
cluster_4_indices = [tf_files_simplified.index(factor) for factor in cluster_4_factors if factor in tf_files_simplified]
cluster_4_matrix = overlap_matrix[:, cluster_4_indices]

# 找到非零行的索引
non_zero_indices = np.where(cluster_4_matrix.any(axis=1))[0]

# 从response_vector1和response_vector2中提取对应的值
selected_values_vector1 = response_vector1[non_zero_indices]
selected_values_vector2 = response_vector2[non_zero_indices]

# 绘制直方图
plot_histogram(selected_values_vector1, 'Histogram of Selected BigWig Signal Vector 1')
plot_histogram(selected_values_vector2, 'Histogram of Selected BigWig Signal Vector 2')


# In[45]:


import matplotlib.pyplot as plt
import numpy as np

def plot_histogram(data, title, bins=50):
    """
    绘制数据的直方图。

    :param data: 数据数组。
    :param title: 直方图的标题。
    :param bins: 直方图的区间数量。
    """
    plt.figure(figsize=(10, 6))
    plt.hist(data, bins=bins, color='blue', alpha=0.7)
    plt.title(title)
    plt.xlabel('Log-transformed value')
    plt.ylabel('Frequency')
    plt.grid(True)
    plt.show()

# 为cluster_4_factors中的每个因子绘制直方图
for factor in cluster_4_factors:
    if factor in tf_files_simplified:
        # 找出因子对应的overlap_matrix中的列
        factor_index = tf_files_simplified.index(factor)
        factor_column = overlap_matrix[:, factor_index]

        # 找到非零行的索引
        non_zero_indices = np.where(factor_column > 0)[0]

        # 从response_vector1和response_vector2中提取对应的值，并计算平均值
        selected_values_vector1 = response_vector1[non_zero_indices]
        selected_values_vector2 = response_vector2[non_zero_indices]
        average_values = (selected_values_vector1 + selected_values_vector2) / 2


        # 绘制经过对数变换的值的直方图
        plot_histogram(average_values, f'Histogram of Log-transformed Average BigWig Signal for {factor}')


# In[46]:


import numpy as np

results = {}  # 用于存储结果的字典

for factor in cluster_4_factors:
    if factor in tf_files_simplified:
        factor_index = tf_files_simplified.index(factor)
        factor_column = overlap_matrix[:, factor_index]
        non_zero_indices = np.where(factor_column > 0)[0]
        selected_values = (response_vector1[non_zero_indices] + response_vector2[non_zero_indices]) / 2
        
        if len(selected_values) > 0:  # 确保数组不为空
            # 使用百分位数作为阈值
            threshold = np.percentile(selected_values, 70)  # 选取高于80%的值（相当于低于20%的阈值）

            # 根据阈值分类
            above_threshold = non_zero_indices[selected_values > threshold]

            # 获取区域的具体位置
            regions_above = [filtered_HMR[i] for i in above_threshold]

            # 存储结果
            results[factor] = {
                'threshold': threshold,
                'above_threshold': regions_above
            }

            print(f"Factor: {factor}")
            print(f"Regions above threshold: {len(regions_above)}")
            print()

# 如果你想查看具体结果，可以直接访问results字典
# 例如，查看某个因子的结果：
for factor, data in results.items():
    print(f"Results for {factor}:")
    print(f"Regions above threshold: {data['above_threshold']}")
    print()


# In[47]:


import numpy as np

results = {}  # 用于存储结果的字典

for factor in cluster_4_factors:
    if factor in tf_files_simplified:
        factor_index = tf_files_simplified.index(factor)
        factor_column = overlap_matrix[:, factor_index]
        non_zero_indices = np.where(factor_column > 0)[0]
        selected_values = (response_vector1[non_zero_indices] + response_vector2[non_zero_indices]) / 2
        
        if len(selected_values) > 0:  # 确保数组不为空
            # 使用百分位数作为阈值
            threshold = np.percentile(selected_values, 30)  # 选取高于80%的值（相当于低于20%的阈值）

            # 根据阈值分类
            above_threshold = non_zero_indices[selected_values < threshold]

            # 获取区域的具体位置
            regions_above = [filtered_HMR[i] for i in above_threshold]

            # 存储结果
            results[factor] = {
                'threshold': threshold,
                'above_threshold': regions_above
            }

            print(f"Factor: {factor}")
            print(f"Regions above threshold: {len(regions_above)}")
            print()

# 如果你想查看具体结果，可以直接访问results字典
# 例如，查看某个因子的结果：
for factor, data in results.items():
    print(f"Results for {factor}:")
    print(f"Regions above threshold: {data['above_threshold']}")
    print()


# In[48]:


import numpy as np

def generate_results_for_factor(factor_name):
    results = {}  # 用于存储结果的字典

    if factor_name in cluster_4_factors and factor_name in tf_files_simplified:
        factor_index = tf_files_simplified.index(factor_name)
        factor_column = overlap_matrix[:, factor_index]
        non_zero_indices = np.where(factor_column > 0)[0]
        selected_values = (response_vector1[non_zero_indices] + response_vector2[non_zero_indices]) / 2

        if len(selected_values) > 0:
            threshold = np.percentile(selected_values, 30)
            above_threshold = non_zero_indices[selected_values < threshold]
            regions_above = [filtered_HMR[i] for i in above_threshold]

            results[factor_name] = {
                'threshold': threshold,
                'above_threshold': regions_above
            }

            print(f"Factor: {factor_name}")
            print(f"Regions above threshold: {len(regions_above)}")
            print()

            return results
    else:
        print("Error: The specified factor is not found in the cluster or file list.")
        return None

# 使用这个函数，例如：
factor_to_check = "MED12_E2_GSE101559_peaks"
result_for_factor = generate_results_for_factor(factor_to_check)
if result_for_factor:
    print(result_for_factor)


# In[49]:


import matplotlib.pyplot as plt

# Values for each cluster
values = [0.261273, 0.033880]
# Corresponding labels for the clusters
labels = ['Cluster 1', 'Cluster 2']

# Create a bar plot
plt.figure(figsize=(8, 6))
plt.bar(labels, values, color=['red', 'green'])

# Adding the value labels on top of each bar
for i, value in enumerate(values):
    plt.text(i, value + 0.01, f'{value:.6f}', ha='center')

# Title and labels
#plt.title('Cluster ESS Values')
plt.ylabel('Average Value')

plt.xticks(labels, labels, fontsize=20)
# Display the plot
plt.tight_layout()
plt.show()


# In[ ]:




