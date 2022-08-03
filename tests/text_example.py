from pybind_sketching import CountMinSketch
import os
import tqdm
import numpy as np
import pandas as pd
from collections import defaultdict # this is just a nice wrapper for a normal dict
import pprint
import matplotlib.pyplot as plt

def main():
    file_name = '../data/shakespeare_complete.txt'
    #file_name = '../data/shakespeare_short.txt'

    num_lines = 0
    word_dict = defaultdict(int)
    with tqdm.tqdm(total=os.path.getsize(file_name)) as progress_bar:
        with open(file_name) as f:
            for line in f:
                progress_bar.update(len(line.encode('utf-8')))
                for word in line.split():
                    num_lines += 1
                    clean_word = ''.join(e for e in word if e.isalnum()) # remove non alphanumeric characters
                    clean_word = clean_word.lower()
                    #print(clean_word)
                    word_dict[clean_word] += 1


    print(f'{num_lines} read in.')
    print(f'{len(word_dict)} words read in.')
    pp = pprint.PrettyPrinter(width=4, compact=True)
    #pp.pprint(word_dict)

    # this bit is cheating but will work for now.
    word_count_list = [None] * len(word_dict)
    word_count_map = {w : {'count' : 0, 'index' : None} for w in word_dict.keys()}
    for idx, (word, count) in enumerate(word_dict.items()):
        word_count_list[idx] = (word, count)
        word_count_map[word]['count'] = count
        word_count_map[word]['index'] = idx


    # Now let's stream the data once more and add the items into the sketch.
    n_buckets = CountMinSketch.suggest_num_buckets(1E-3)
    n_hashes = CountMinSketch.suggest_num_hashes(1-1E-4)
    sk = CountMinSketch(n_hashes, n_buckets, 100)


    with tqdm.tqdm(total=os.path.getsize(file_name)) as progress_bar:
        with open(file_name) as f:
            for line in f:
                progress_bar.update(len(line.encode('utf-8')))
                for word in line.split():
                    num_lines += 1
                    clean_word = ''.join(e for e in word if e.isalnum())  # remove non alphanumeric characters
                    clean_word = clean_word.lower()
                    clean_word_loc = word_count_map[clean_word]['index']
                    sk.update(clean_word_loc)

    for word in word_count_map.keys():
        word_count_map[word]['estimate'] = sk.get_estimate(word_count_map[word]['index'])
    #pp.pprint(word_count_map)
    shakespeare_df = pd.DataFrame.from_dict(word_count_map, orient='index')
    shakespeare_df.sort_values(by="count", ascending=False, inplace=True)
    shakespeare_df.reset_index(inplace=True)
    shakespeare_df.rename(columns={'level_0' : 'word', 'index' : 'map_index'}, inplace=True)
    shakespeare_df.reset_index(inplace=True)
    print(shakespeare_df.head(15))
    print(sk)
    table = np.array(sk.get_table())
    np.save("cm_sketch.npy", table)
    fig, ax = plt.subplots()
    shakespeare_df.plot.scatter(x='index', y="count", ax=ax, c="blue", s=1.0, marker='.', label='Count')
    shakespeare_df.plot.scatter(x='index', y="estimate", ax=ax, c="red", s=1.0, marker='^', label='Estimate')
    ax.legend()
    ax.grid()
    ax.set_yscale('log', base=10)
    ax.set_xscale('log', base=10)
    plt.show()


    # print("{:>6} {:8} {:>5}".format('Hashes', 'Buckets', 'Mass'))
    # print("{:>6} {:8} {:>5}".format(sk.get_num_hashes(), sk.get_num_buckets(), sk.get_total_weight()))


if __name__ == '__main__':
    main()