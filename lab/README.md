# Interactive lab

The browser lab is a self-contained static application. Open `index.html` via
GitHub Pages or a local HTTP server; its React and case-data assets are stored
in this same directory.

The source uses in-browser Babel compilation for this initial public release.
The numerical reference package and reproducibility documentation are in
`../code/`.

Use the `中文 / English` switch in the upper-right corner. The selected
language is retained in local storage and can also be set with `?lang=en` or
`?lang=zh`.

Application Case B is event-triggered: the slow vehicle and fixed upstream VSL
both start at `t = 300 s`. The VSL is strictly off before the event and ramps in
over 30 seconds.
