# Contributing

We welcome contributions!

## Development Setup

1.  Clone the repository:
    ```bash
    git clone https://github.com/elkins/synth-pdb.git
    ```
2.  Install in editable mode with dev dependencies:
    ```bash
    pip install -e ".[dev]"
    ```
3.  Run tests:
    ```bash
    pytest
    ```

## Code Style

We use `black` and `ruff` for linting.

```bash
ruff check .
black .
```
