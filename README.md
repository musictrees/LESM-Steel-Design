# LESM - Linear Elements Structure Model

<p align=center><img height="100.0%" width="100.0%" src="https://gitlab.com/rafaelrangel/lesm/-/raw/master/docs/images/logos/logo_lesm_landscape.png"></p>

LESM is a graphical-interactive program for structural analysis of models composed of linear elements,
i.e. uniaxial elements with one dimension much larger than the others,
such as bars (axial behavior only) and beams (axial, flexural, and torsion behaviors).

Its purpose is to offer a free and user-friendly tool to be used both in academia and industry.
It is intended to help students, teachers, engineers and architects improve the quality of their work.

For more information, visit its [website][lesm_website] and [YouTube channel][youtube_channel].

If you have used the program, pelase leave your [feedback][feedback_link].

## Table of Contents
- [Main Features](#main-features)
- [Implementation Aspects](#implementation-aspects)
- [Instructions](#instructions)
    - [Running the Source Code in MATLAB](#running-the-source-code-in-matlab)
    - [Testing](#testing)
	- [Installing the MATLAB App](#installing-the-matlab-app)
	- [Installing the Standalone Executable](#installing-the-standalone-executable)
- [Examples](#examples)
- [Documentation](#documentation)
- [History](#history)
- [How to Contribute](#how-to-contribute)
- [How to Cite](#how-to-cite)
- [Authorship](#authorship)
- [License](#license)

## Main Features

The available **model types** are:

- Two-dimensional (plane) truss
- Two-dimensional (plane) frame
- Three-dimensional (spatial) truss
- Three-dimensional (spatial) frame
- Grillage

The available **analysis types** are:

- Linear-elastic static
- Linear-elastic dynamic (time-domain)

The program has an easy-to-use **interactive GUI** (Graphical User Interface) for modeling, pre-processing, and post-processing the analysis results.

It can be used by running its source code in the **MATLAB environment**, or by installing the **MATLAB App** or the **standalone executable**.

## Implementation Aspects

LESM is entirely developed in [MATLAB][matlab_website].

The analysis solver of the program is implemented using the Object-Oriented Programming (OOP) paradigm to offer modularity and extensibility.

The GUI is developed in the Graphical User Interface Development Environment (GUIDE) of MATLAB.

## Instructions

LESM is used via a graphical interface, which can be launched by directly running its source code in MATLAB, or by using the MATLAB App or the standalone executable after installing them.

Tutorials on how to use the interface can be found in the [manuals][manuals_folder_link], on the [Wiki page][wiki_link], and in the [YouTube channel][youtube_channel].

**Attention**:
The minimum recommended screen resolution for using the interface is 1366x768, and preference should be given to a resolution of 1920x1080.
A lower resolution can lead to clipping or misalignment of some interface components.

It is also recommended to use the interface maximized.

### Running the Source Code in MATLAB

To run LESM in any MATLAB version, execute the script file [*lesm.m*][lesm_file_link] located inside the folder [*src*][src_folder_link] to launch the interface.

### Testing

When using LESM in MATLAB, recursive tests are available to verify that the solver of the program is working correctly and that the current results are matching with the reference results.

To **run the tests**, execute the script file [*run_tests.m*][run_tests_link] located inside the folder [*tests*][tests_folder_link].
In that file, you can set the tolerance for comparing the current results with the reference results. A dialog box will pop up to select the input files of the tests to be run, which are located inside the sub-folder [*test_models*][test_models_link].
The result of each test is then printed in the MATLAB command window.

To **generate or update the reference results** of the test models, execute the script file [*update_results.m*][upd_tests_link] located inside the folder [*tests*][tests_folder_link].
In that file, you can set the number of decimal places to print the reference results A dialog box will pop up to select the input files of the tests to be updated, which are located inside the sub-folder [*test_models*][test_models_link].
It will run the selected tests and overwrite existing reference results.

### Installing the MATLAB App

This option adds LESM as a new application in the MATLAB environment.

The MATLAB App installer file can be found in the folder [installation][install_folder_link] (file *.mlappinstall*).
To install the MATLAB App of LESM, double-click the installer file from the MATLAB Current Folder browser.

### Installing the Standalone Executable

**Available only for Windows and MacOS systems.**

This option is indicated for those who do not have MATLAB installed.

It will install the MATLAB runtime and the executable standalone of LESM.
The runtime is a set of libraries that make it possible to run compiled applications without the need of having MATLAB installed or the license to install it.
If you already have MATLAB on your machine, the installer will not install the runtime, only the executable standalone of LESM.
Otherwise, the installer will download the MATLAB runtime, so **internet connection** may be required.

The files of the LESM installer can be found in the folder [installation][install_folder_link].
To start the installation process, run the appropriate installer file for your operating system (*.exe* or *.app*).

The installation process is straightforward.
It is described [here][install_guide_link] and illustrated in the images below for a Windows system.

<p float="left">
<img src="https://gitlab.com/rafaelrangel/lesm/-/raw/master/docs/installation/win_install_01.png" width="275"/>
&nbsp;
<img src="https://gitlab.com/rafaelrangel/lesm/-/raw/master/docs/installation/win_install_02.png" width="275"/>
&nbsp;
<img src="https://gitlab.com/rafaelrangel/lesm/-/raw/master/docs/installation/win_install_03.png" width="275"/>
</p>

<p float="left">
<img src="https://gitlab.com/rafaelrangel/lesm/-/raw/master/docs/installation/win_install_04.png" width="275"/>
&nbsp;
<img src="https://gitlab.com/rafaelrangel/lesm/-/raw/master/docs/installation/win_install_05.png" width="275"/>
&nbsp;
<img src="https://gitlab.com/rafaelrangel/lesm/-/raw/master/docs/installation/win_install_06.png" width="275"/>
</p>

## Examples

A collection of sample models of different types is available inside the folder [*examples*][examples_folder_link].

Each model is stored in a file with the extension *.lsm*, which is saved and loaded by the program.

The sample models can be loaded through the graphical interface of LESM to explore the features of the program.

<p align=center>
<img height="100.0%" width="35.0%" src="https://gitlab.com/rafaelrangel/lesm/-/raw/master/docs/images/examples/example_01.png">
&nbsp;
<img height="100.0%" width="30.0%" src="https://gitlab.com/rafaelrangel/lesm/-/raw/master/docs/images/examples/example_02.png">
&nbsp;
<img height="100.0%" width="30.0%" src="https://gitlab.com/rafaelrangel/lesm/-/raw/master/docs/images/examples/example_03.png">
</p>

<p align=center>
<img height="100.0%" width="35.0%" src="https://gitlab.com/rafaelrangel/lesm/-/raw/master/docs/images/examples/example_04.png">
&nbsp;
<img height="100.0%" width="25.0%" src="https://gitlab.com/rafaelrangel/lesm/-/raw/master/docs/images/examples/example_05.png">
&nbsp;
<img height="100.0%" width="35.0%" src="https://gitlab.com/rafaelrangel/lesm/-/raw/master/docs/images/examples/example_06.png">
</p>

## Documentation

* **Book**: The book [*Analise Matricial de Estruturas com Orientacao a Objetos*][book_link] (Matrix Structural Analysis with Object-Oriented Programming) describes the theory and implementation
of the direct stiffness method using LESM as a supporting program.
The source code of the LESM version that accompanies the book is fully documented through [HTML files] and [UML diagrams].

* **Wiki pages**: General and important information for starting using the program is available in the [Wiki pages][wiki_link].

* **Manuals**: Manuals explaining the considerations of the program and the use of the graphical interface for different release versions are available [here][manuals_folder_link].

* **YouTube**: Video-tutorials and presentations in congresses and conferences can be watched in the [YouTube channel][youtube_channel]

## History

- **LESM 3.0 (Dynamic version)**: May 2022
  
  Developed during an internship project by Pedro Lopes,
  advised by Luiz Fernando Martha and Rafael Rangel -
  Tecgraf/PUC-Rio.
  
  Consolidated during the PhD thesis of Claudio Horta, advised by Luiz Fernando Martha -
  Department of Civil and Environmental Engineering, PUC-Rio, 2022.
  
  Main new features:
  Dynamic analysis option.

- **LESM 2.0 (Graphical version)**: October 2018
  
  Developed during an internship project by Pedro Lopes,
  advised by Luiz Fernando Martha and Rafael Rangel -
  Tecgraf/PUC-Rio.
  
  Main new features:
  Development of a graphical interface with mouse interaction for modelling and post-processing,
  inclusion of load cases and combinations, inclined and spring supports, and semi-rigid connections.

- **LESM 1.0 (Book version)**: May 2017
  
  Initial version documented in the book
  [*Analise Matricial de Estruturas com Orientacao a Objetos*][book_link] (Matrix Structural Analysis with Object-Oriented Programming),
  by Luiz Fernando Martha.
  
  Developed in the undergraduate thesis
  [Development of a Graphic Program for Structural Analysis of Linear Element Models][thesis_rafael_link],
  by Rafael Rangel, advised by Luiz Fernando Martha -
  Department of Civil and Environmental Engineering, PUC-Rio, 2016.
  
  Extended in the undergraduate thesis
  [Extensao de Programa Grafico para Analise de Trelicas e Porticos Espaciais via MATLAB e GUI][thesis_murilo_link],
  by Murilo Felix, advised by Luiz Fernando Martha -
  Department of Civil and Environmental Engineering, PUC-Rio, 2017.

## How to Contribute

Please check the [contribution guidelines][contribute_link].

## How to Cite

When citing LESM in your work, please refer to one of the technical journal publications listed on its [website][publications_link].

## Authorship

- **Luiz Fernando Martha** (<lfm@tecgraf.puc-rio.br>)
- **Rafael L. Rangel** (<rafaelrangel@tecgraf.puc-rio.br>)
- **Pedro C. Lopes** (<cortezpedro@tecgraf.puc-rio.br>)
- **Cláudio Horta** (<claudiohorta@tecgraf.puc-rio.br>)

Pontifical Catholic University of Rio de Janeiro (PUC-Rio) -
[Department of Civil and Environmental Engineering][civil_website]

Tecgraf Institute of Technical-Scientific Software Development of PUC-Rio
([Tecgraf/PUC-Rio][tecgraf_website])

<p float="left">
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
<img src="https://gitlab.com/rafaelrangel/lesm/-/raw/master/docs/images/logos/logo_puc.png" width="300"/>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
<img src="https://gitlab.com/rafaelrangel/lesm/-/raw/master/docs/images/logos/logo_tecgraf.png" width="400"/> 
</p>

## License

LESM is licensed under the [MIT license][mit_license_link],
which allows the program to be freely used by anyone for modification, private use, commercial use, and distribution, only requiring preservation of copyright and license notices.

No liability and warranty is provided.
Neither the authors, nor PUC-Rio, nor Tecgraf/PUC-Rio, nor any related institution are responsible for any use or misuse of the program and the results.
The aforementioned assume no liability or responsibility to any person or company for direct or indirect damages resulting from the use of any information or the use of any of the software made available here.
The user is responsible for any and all conclusions made while using the program.

[lesm_website]:         https://web.tecgraf.puc-rio.br/lesm
[youtube_channel]:      https://www.youtube.com/channel/UCoxRpAftcAhfqREELS0kj6w
[feedback_link]:        https://forms.gle/tyq7VXX2n2EGwFST9
[matlab_website]:       https://www.mathworks.com/
[manuals_folder_link]:  https://gitlab.com/rafaelrangel/lesm/-/tree/master/docs/manuals
[wiki_link]:            https://gitlab.com/rafaelrangel/lesm/-/wikis/home
[lesm_file_link]:       https://gitlab.com/rafaelrangel/lesm/-/blob/master/src/lesm.m
[src_folder_link]:      https://gitlab.com/rafaelrangel/lesm/-/tree/master/src
[run_tests_link]:       https://gitlab.com/rafaelrangel/lesm/-/blob/master/tests/run_tests.m
[tests_folder_link]:    https://gitlab.com/rafaelrangel/lesm/-/tree/master/tests
[test_models_link]:     https://gitlab.com/rafaelrangel/lesm/-/tree/master/tests/test_models
[upd_tests_link]:       https://gitlab.com/rafaelrangel/lesm/-/blob/master/tests/update_results.m
[install_folder_link]:  https://gitlab.com/rafaelrangel/lesm/-/tree/master/docs/installation
[install_guide_link]:   https://www.mathworks.com/help/compiler/install-deployed-application.html
[examples_folder_link]: https://gitlab.com/rafaelrangel/lesm/-/tree/master/examples
[book_link]:            https://www.grupogen.com.br/e-book-analise-matricial-de-estruturas-com-orientacao-a-objetos
[HTML files]:           https://web.tecgraf.puc-rio.br/lesm/v1/publish/main.html
[UML diagrams]:         https://web.tecgraf.puc-rio.br/lesm/v1/uml/Astah_html/index.html
[thesis_rafael_link]:   https://web.tecgraf.puc-rio.br/~lfm/teses/RafaelRangel-MonografiaGraduacao-2016.pdf
[thesis_murilo_link]:   https://web.tecgraf.puc-rio.br/~lfm/teses/MuriloFelixFilho-MonografiaGraduacao-2016.pdf
[contribute_link]:      https://gitlab.com/rafaelrangel/lesm/-/blob/master/CONTRIBUTING.md
[publications_link]:    https://web.tecgraf.puc-rio.br/lesm/publications.html
[civil_website]:        https://bananastudio.com.br/civil/web/
[tecgraf_website]:      https://www.tecgraf.puc-rio.br/
[mit_license_link]:     https://choosealicense.com/licenses/mit/
