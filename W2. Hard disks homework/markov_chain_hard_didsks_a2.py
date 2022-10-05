"""
The solution of the A2 part:
reimplmentation of Markov chain sampling
of 4 hard disks deposition in the box

required python 3.9
"""

from copy import deepcopy
import random
import math

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


def markov_disks_box_iterator(initial_state, sigma=0.15, delta=0.01):
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


def run_multi_markov(initial_state, n_steps=1_000_000, sigma=0.15, step=0.1, del_xy=0.05):
    configurations = {
        "conf_a": A,
        "conf_b": B,
        "conf_c": C,
    }
    hits = {"conf_a": 0, "conf_b": 0, "conf_c": 0}
    samples = []
    markov_disks_box_generator = markov_disks_box_iterator(
        initial_state,
        sigma,
        step,
    )
    
    for _ in range(n_steps):
        x_vec = next(markov_disks_box_generator)
        samples.append(x_vec)
        for conf in configurations:
            condition_hit = True
            for b in configurations[conf]:
                condition_b = min(
                    max(
                        abs(a[0] - b[0]),
                        abs(a[1] - b[1])
                    ) for a in x_vec
                ) < del_xy
                condition_hit *= condition_b
            if condition_hit:
                hits[conf] += 1

    for conf in configurations:
        print(conf, hits[conf])
    return samples

init_state = [[0.25, 0.25], [0.75, 0.25], [0.25, 0.75], [0.75, 0.75]]

for nsteps in [10_000, 100_000, 1_000_000, 10_000_000]:
    print(f"{nsteps} steps:")
    for run in range(3):
        print(f"Run #{run + 1}")
        run_multi_markov(
            init_state,
            nsteps,
            sigma=0.15,
            step=0.1,
            del_xy=0.05,
        )
        print()