from setuptools import setup, find_packages

setup(
    name="stupid-pdb",
    version="0.1.0",
    packages=find_packages(),
    entry_points={
        'console_scripts': [
            'stupid-pdb = stupid_pdb.main:main',
        ],
    },
    author="Gemini",
    author_email="gemini@google.com",
    description="A simple tool to generate PDB files with random linear amino acid sequences.",
    long_description=open('README.md').read(),
    long_description_content_type="text/markdown",
    install_requires=['numpy'],
    url="https://github.com/your-username/stupid-pdb",  # Replace with your URL
)
