# CHORD Protocol in MPI

This project implements a simplified version of the **CHORD distributed hash table (DHT)** protocol using MPI. Each MPI process represents a node in a static CHORD ring and participates in distributed key lookups using **finger tables** for logarithmic routing.

## Overview

- Each node has a unique **CHORD ID**, a **successor**, a **predecessor**, and a **finger table**.
- Lookup requests for keys are routed through the network using `closest_preceding_finger()` to achieve **O(log N)** routing.
- Communication is handled with `MPI_Send` and `MPI_Recv`.
- The system is **static**: nodes do not join or leave during execution.

## Features

- Build finger tables (`build_finger_table()`) for efficient routing.
- Logarithmic routing using `closest_preceding_finger()`.
- Distributed lookup handling (`handle_lookup_request()`).
- Circular ID space with consistent hashing.
- Service loop for handling requests, responses, and termination (`TAG_DONE`).

## Technologies and Concepts

- **MPI (Message Passing Interface):** For inter-process communication between CHORD nodes.  
- **Distributed Hash Tables (DHTs):** Storing and retrieving keys across a distributed network.  
- **Finger Tables:** Optimized routing table to achieve logarithmic lookups.  


## How to run

The checker verifies correctness and efficiency.

```bash
cd checker/
./checker.sh
```

