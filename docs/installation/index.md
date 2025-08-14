---
title: Installing PySTARS
nav_order: 2
---

# Installing PySTARS

This guide walks you through setting up PySTARS on your system.

---

## 1. Clone the Repository

```bash
git clone https://github.com/POptUS/RanDFO.git
cd RanDFO/PySTARS
```

---

## 2. Create a Virtual Environment

It’s recommended to use a virtual environment to keep dependencies isolated.

```bash
python -m venv venv
```

### Activate the Virtual Environment

- **Windows (Command Prompt):**
  ```cmd
  venv\Scripts\activate
  ```

- **Windows (Git Bash / PowerShell):**
  ```bash
  source venv/Scripts/activate
  ```

- **Linux / macOS:**
  ```bash
  source venv/bin/activate
  ```

---

## 3. Install Required Dependencies

PySTARS requires the following Python packages:

| Package      | Purpose                               |
|--------------|---------------------------------------|
| `setuptools` | Packaging and installation support    |
| `pandas`     | Data handling and logging             |
| `matplotlib` | Plotting convergence graphs           |
| `scipy`      | Scientific computations               |
| `joblib`     | Parallel processing                   |

You can install them all using the included **requirements.txt**:

```bash
pip install -r requirements.txt
```

[View requirements.txt](https://github.com/POptUS/RanDFO/tree/add-pystars/PySTARS/requirements.txt)
