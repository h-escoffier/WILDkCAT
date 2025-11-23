# Installation

## Setup a Python environment 

To ensure a clean and isolated setup, we recommend to use [uv](https://docs.astral.sh/uv/), a lightweight tool that simplifies Python environment and package management. If you don’t have it yet:

=== "macOS / Linux"
    ```bash
    curl -LsSf https://astral.sh/uv/install.sh | sh
    ```

=== "Windows"
    ```powershell
    powershell -ExecutionPolicy ByPass -c "irm https://astral.sh/uv/install.ps1 | iex"
    $env:Path += ";$env:USERPROFILE\.local\bin"
    ```

Create and activate a virtual environment with uv:

=== "macOS / Linux"
    ```bash
    uv venv
    source .venv/bin/activate
    ```

=== "Windows"
    ```powershell
    uv venv
    .venv\Scripts\activate
    ```

!!! tip 

    Using **uv** helps manage dependencies easily and keeps your environment clean, avoiding conflicts with other Python packages.

---

## Install WILDkCAT from PyPI 

Install [WILDkCAT](https://pypi.org/project/wildkcat/) directly from PyPI:

```bash
uv pip install wildkcat 
```

Verify the installation by checking the version:

```bash 
wildkcat --version
```

---

### Environment Setup

Provide your **BRENDA login credentials** and **Entrez API email adress** to query the BRENDA enzyme database and NCBI database.

Create a file named `.env` in the root of your project with the following content:

```bash
ENTREZ_EMAIL=your_email@example.com
BRENDA_EMAIL=your_email@example.com
BRENDA_PASSWORD=your_password
```

!!! note

    - Replace the placeholders with the credentials from the [BRENDA website](https://www.brenda-enzymes.org) account you created.

!!! danger 

    - Keep the file `.env` **private** (e.g., add it to your `.gitignore`) as it contains sensitive information.

--- 

## Install CataPro

[CataPro](https://github.com/zchwang/CataPro?tab=readme-ov-file#create-the-catapro-environment) is required for predicting kcat values using machine learning

### 1. Install dependencies

```bash
uv pip install torch transformers numpy pandas RDKit sentencepiece
```

### 2. Clone the CataPro repository

```bash
git clone https://github.com/zchwang/CataPro
```

### 3. Set up Git LFS 

CataPro uses [Git Large File Storage (LFS)](https://git-lfs.github.com/) to handle large model files. 
If you don't have Git LFS installed, you can install it using the following command:

```bash
git lfs install
```

### 4. Download the models

CataPro requires the models **ProtT5 model** and **MolT5 model**.

!!! warning 

    The models prot_t5_xl_uniref50 and molt5-base-smiles2caption required for CataPro are 64 and 1.9 GB, respectively. 


=== "macOS / Linux"
    ```bash
    cd CataPro/models/

    LFS_SKIP_SMUDGE=1 git clone https://huggingface.co/Rostlab/prot_t5_xl_uniref50
    cd prot_t5_xl_uniref50
    git lfs pull

    cd ..

    LFS_SKIP_SMUDGE=1 git clone https://huggingface.co/laituan245/molt5-base-smiles2caption
    cd molt5-base-smiles2caption
    git lfs pull

    cd ../..
    ```

=== "Windows"
    ```powershell
    git -c filter.lfs.smudge= -c filter.lfs.required=false clone https://huggingface.co/Rostlab/prot_t5_xl_uniref50
    cd prot_t5_xl_uniref50
    git lfs pull

    cd ..

    git -c filter.lfs.smudge= -c filter.lfs.required=false clone https://huggingface.co/laituan245/molt5-base-smiles2caption
    cd molt5-base-smiles2caption
    git lfs pull

    cd ../..
    ```


---
