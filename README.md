# PV-MSA
A pure JavaScript wrapper combining BioJS MSA and PV.

### Getting started:

To use, import the required scripts and CS into your HTML page:

```
<head>
    <link href="/lib/msa/msa.css" rel="stylesheet">
    <link href="/lib/pv-msa/pv-msa.css" rel="stylesheet">
</head>

...

<script src="/lib/pv/bio-pv.min.js"></script>
<script src="/lib/msa/msa.min.js"></script>
<script src="/lib/pv-msa/pv-msa.js"></script>
```

Somewhere in your HTML file, create an empty `div` that will be used to house the protein viewer:

```
<div id="viewer"></div>
```

In you JavaScript, create the protein viewer and load the structure and sequence into it:

```javascript
// create the viewer
var viewer = new PV("viewer");
viewer.render();

// create a structure object
var seq1 = new Structure(
    "target001", // id
    "target001", // name
    "VLSPADKTNVKAAWGKVGAHAG--EYGAEALERMFLSFPTTKTYFP---HFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR", // sequence
    "/lib/data/target001.pdb", // URL for PDB file
    "red", // color
    1,
    false
)

// add structure to viewer
viewer.addStructure(seq1, true);
```

### Demo

*Note: you will need to have npm and node installed for the demo to work.*

An example of how to use PV-MSA can be found in `index.html`. To run a demo, install express using the following command:

```
npm install
```


Once express is installed, you can easily run the demo by executing:

```
npm start
```
