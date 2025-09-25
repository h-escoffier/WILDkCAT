# Installation

## Install WILDkCAT

### Install from PyPI

Install [WILDkCAT](https://pypi.org/project/wildkcat/) directly from PyPI:

```bash
pip install wildkcat
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

To predict missing kcat values using machine learning, install [CataPro](https://github.com/zchwang/CataPro?tab=readme-ov-file#create-the-catapro-environment).

CataPro also requires two additional pre-trained models: [prot_t5_xl_uniref50](https://huggingface.co/Rostlab/prot_t5_xl_uniref50) and [molt5-base-smiles2caption](https://huggingface.co/laituan245/molt5-base-smiles2caption).
