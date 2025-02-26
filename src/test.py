from code_snippet_FC import sample_from_mixture, KSD_loss_RBF
import re
import numpy as np

n = 10000

def subset_from_log(filename):
    ans = None
    pattern = re.compile(r'^\s*Active points:\s*(.*)$')
    
    with open(filename, 'r', encoding='utf-8') as file:
        for line in file:
            match = pattern.match(line)
            if match:
                ans = list(map(int, match.group(1).split()))
    return ans


points = sample_from_mixture(n, 2, seed=42)

subset = subset_from_log('output10000/460.txt')

sub = points[subset]
v = KSD_loss_RBF(sub, 1, len(subset), 2)
vv = KSD_loss_RBF(points, 1, n, 2)
subb = vv[0].numpy().astype(np.float64)[np.ix_(subset, subset)]
print(subb.mean())
print(v.mean())

