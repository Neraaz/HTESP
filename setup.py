"""Setup script for installation of HTESP package"""
import glob
import setuptools
setuptools.setup(
    name="HTESP",
    version="v1.0",
    author="Niraj K. Nepal",
    author_email="nepalneeraz@gmail.com",
    license="",
    description="Package to perform high-throughput QE and VASP calculations",
    url="https://github.com/Neraaz/HTESP",
    install_requires=[
        "scipy",
        "numpy",
        "pandas",
        "pymatgen",
        "mp_api",
        "ase",
        "bsym",
        "lmfit",
        "qmpy-rester",
        "ifermi",
    ],
    # scripts=["alignn/alignn_train_folder.py"],
    packages=setuptools.find_packages(),
    python_requires=">=3.10",
    classifiers=[
        "Programming Language :: Python :: 3 + Bash Scripting",
        "Operating System :: UNIX",
        "License :: OSI Approved :: MIT License",
    ],
    scripts=[f for f in glob.glob('src/*.py') if not f.endswith('__init__.py')],
    keywords="Superconductivity,QE,VASP,WANNIER90,WannierTools,Phonopy,EPW,el-ph coupling",
    entry_points={
        'console_scripts': ['mainprogram = src.mainprogram:main',],
    }
)
