"""
The solution of the B2 part:
reimplmentation of Markov chain sampling
of 4 hard disks deposition in the box

required python 3.9
"""

from copy import deepcopy
import random
from sys import argv
import matplotlib.pyplot as plt
import numpy as np


# Lookup satates
A = (
    (0.30, 0.30),
    (0.30, 0.70),
    (0.70, 0.30),
    (0.70, 0.70),
)
B = (
    (0.20, 0.20),
    (0.20, 0.80),
    (0.75, 0.25),
    (0.75, 0.75),
)

C = (
    (0.30, 0.20),
    (0.30, 0.80),
    (0.70, 0.20),
    (0.70, 0.70),
)


def markov_disks_box_iterator(initial_state, sigma=0.15, delta=0.1):
    L = deepcopy(initial_state)
    yield list(L)
    while True:
        a_index = random.randint(0, 3)
        a = L[a_index]
        b = [
            a[0] + random.uniform(-delta, delta),
            a[1] + random.uniform(-delta, delta),
        ]
        min_dist = min(
            pow(b[0] - L[c][0], 2) + pow(b[1] - L[c][1], 2)
            for c in range(4) if c != a_index
        )
        box_cond = min(b[0], b[1]) < sigma or max(b[0], b[1]) > 1.0 - sigma
        if not (box_cond or min_dist < 4.0 * sigma ** 2):
            L[a_index] = b
        yield list(L)  # list(L) freezes the current state


init_state = [[0.25, 0.25], [0.75, 0.25], [0.25, 0.75], [0.75, 0.75]]
n_steps = int(argv[1])

samples = np.asarray([
    l
    for _, l in zip(
        range(n_steps),
        markov_disks_box_iterator(init_state, sigma=0.1197)
    )
])

fig, ax = plt.subplots(figsize=(8,7))
ax.hist(
    samples[:, 0, 0],
    bins=np.linspace(0, 1, 120),
    color="k",
)
ax.set_xlabel('x')
ax.set_ylabel('frequency')
ax.set_title('Markov sampling: x coordinate histogram (density eta=0.18)')
fig.savefig('markov_disks_histo.png')
