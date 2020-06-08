# ciuGUI—a simple tool for comparing two CIU datasets

Two .csv files with the standard CIU format are imported--these must contain the same number of rows and columns. The data can then be smoothed with either 1D area or intensity normalisation.

The data are displayed as a four-pane figure. The top two panes show the heat-maps for each data set in the selected colormap. The bottom-left pane shows the RMSD heat-map and the rolling RMSD (across the collision energies). The bottom-right pane shows the overlaid mobiligrams for a selected energy, the corresponding energy is also displayed as a line on the top panes.

This is a free and open-source software for the specific task of comparing CIU datasets. Other, more comprehensive, software is available such as [CIUSuite](https://github.com/RuotoloLab/CIUSuite2), [ORIGAMI](https://github.com/lukasz-migas/ORIGAMI) and [PULSAR](http://pulsar.chem.ox.ac.uk/).

```
Copyright (c) 2016-2020 Jedd Bellamy-Carter, University of Nottingham
Released under a GNU Public License
```

## Usage

Simply run from the command line:
```bash
python ciuGui.py
```
If you prefer, create a desktop shortcut to run the above command.

![Example Screenshot](https://github.com/jbellamycarter/ciuGUI/blob/master/example.jpg)

### Keyboard Shortcuts
|              |                                 |
|--------------|---------------------------------|
|`<Control>+O` | Open CIU Data files             |
|`<Return>`    | Update Graphs                   |
|`<Control>+S` | Save Processed CIU Data to file |
|`<Control>+E` | Export Image to file            |

### Step-by-step
1. Once launched, you can select `Open CIU Data` to choose two CIU data sets, these **must** have the same dimensions otherwise they will not be plotted.
2. You can then give each dataset a name and list the charge state—this software uses $E_{lab}$ values.
3. (Optional) You can change the smoothing parameters and colormap (following `matplotlib` names) in use.
4. Select a collision energy to examine and `Update Graphs`—`<Return>` keyboard shortcut.


#### Thank You for using ciuGUI!
