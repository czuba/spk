# spk

Ephys spike loading &amp; analysis tools (Matlab, PLDAPS, Plexon, Kilosort...)

- Analysis functions currently coded up specific to mtSpeedVar project, but will eventually become a more generalized repo.



### Dependencies

- Loading .PDS data files requires a current version of the **[PLDAPS (glDraw branch)](https://github.com/czuba/PLDAPS.git)** in the Matlab path.
_(...without this, custom data classes, like the `condMatrix`, cannot be loaded properly)_


- All you should need in your path for these analyses to work is PLDAPS and this SPK repo