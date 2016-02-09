var msa = require("msa");

var temp;
function ResidueRange(start, end) {
    this.start = start;
    this.end = end;
}

function Residue(name, pos) {
    this.name = name;
    this.pos = pos;
}

function Structure(id, name, seq, struc_url, color, height, ref) {
    this.id = id;
    this.name = name;
    this.seq = seq;
    this.struc_url = struc_url;
    this.color = color;
    this.height = height;
    this.ref = ref;
    
    this.residues = [];
    
    var count = 0;
    for(var i = 0; i < this.seq.length; i++) {
        var name = this.seq[i];
        var pos = -1;
        
        if(name != 'X' && name != '-') {
            count++;
            pos = count;
        }
            
        this.residues.push(new Residue(name, pos));
    }
    
    this.struc = null;
    this.geom = null;
    this.seq_obj = null;
}

function PV(viewer_id) {
    var self = this;
    
    self.parent = null;
    if (viewer_id != null) {
        self.parent = document.getElementById(viewer_id);
    }
    
    self.viewer_container = null;
    self.viewer = null;
    self.geom = null;
    self.prevPicked = null;
    self.reference = null;
    
    self.atom_name = null;
    
    self.alignment_container = null;
    self.alignment = null;
    
    self.mask = null;
    
    self.options = {
        width: 800,
        seq_height: 100,
        struc_height: 600,
        struc_quality: "medium"
    };
    
    /*self.options = {
        width: 800,
        alignment: {
            height: 80,
            colorBy: 'default'
        },
        structure: {
            height:600,
            quality: "medium"
        }
    }*/
    
    self.structures = [];
    self.columns = null;
        
    self.render = function() {
        /* clear parent div of current contents */
        self.parent.innerHTML = "";
        self.parent.style.width = self.options.width + 1;
        
        /* create the div for PV and MSA*/
        self.viewer_container = document.createElement("div");
	    self.viewer_container.id = "pv_" + self.parent.id;
        self.viewer_container.style.width = self.options.width + 1;
	    
	    self.atom_name = document.createElement("div");
	    self.atom_name.id = "picked-atom-name";

        self.alignment_container = document.createElement("div");
	    self.alignment_container.id = "msa_" + self.parent.id;
        
        self.parent.appendChild(self.viewer_container);
        self.parent.appendChild(self.atom_name);
        self.parent.appendChild(self.alignment_container);
        
        self.loadMSA();
        self.loadViewer();
    }
    
    self.renderSeparate = function(viewer, msa, atom_name) {
        self.parent = null;
        self.viewer_container = document.getElementById(viewer);
        self.viewer_container.style.width = self.options.width + 1;
        
        self.atom_name = document.getElementById(atom_name);
        
        self.alignment_container = document.getElementById(msa);
        self.loadMSA();
    }
    
    self.setViewerHeight = function(height) {
        self.viewer_container.style.height = height;
        self.options.struc_height = height;
        self.viewer.resize(self.options.width-3, height);
    }
    
    self.setWidth = function(width) {
        self.options.width = width;
        
        self.parent.style.width = width;
        self.viewer_container.style.width = width + 1;
        
        self.viewer.resize(width-3, self.options.struc_height);
        self.alignment.g.zoomer.set("alignmentWidth", width - 100);
    }
    
    self.addStructure = function(structure, show) {
        self.structures.push(structure);
        
        if(show) {
	   	    self.loadViewer();
            self.showStructure(structure);
            self.showSequence(structure);
        }
        
        if(self.columns == null)
            self.columns = new Array(structure.residues.length);
            
        return structure;
    }
    
    self.showSequence = function(structure) {
        if(structure.seq_obj != null) {
            structure.seq_obj.set("hidden", false);
        } else {
            self.alignment.seqs.add(structure);
            structure.seq_obj = self.alignment.seqs.get(structure.id);
        }
    }
    
    self.hideSequence = function(structure) {
        if(structure.seq_obj != null) {
            structure.seq_obj.set("hidden", true);
        }
    }
    
    self.showStructure = function(structure) {
        self.showLoading();
        
        if(structure.struc == null) {
    		pv.io.fetchPdb(structure.struc_url, function(s) {
    			structure.struc = s;
    		    
    	   	    self.showStructure(structure);
          	});
    	} else {
	   	    if(self.reference != null) {
                out = pv.mol.superpose(structure.struc, self.reference);
            } else {
                self.reference = structure.struc;
            }
            
            if(structure.geom == null) {  
        	    structure.geom = self.viewer.cartoon(structure.name, structure.struc, { color : color.uniform(structure.color) });
        		self.viewer.centerOn(structure.struc);
        		
                self.viewer.autoZoom();
                
                for(var i = 0; i < self.columns.length; i++) {
                    if(self.columns[i] == true)
                        self.setResidueSelected(structure, i, true);
                }
            } else {
                self.viewer.show(structure.name);
                self.viewer.requestRedraw();
            }
            
            self.hideLoading();
    	}
    }
    
    self.hideStructure = function(structure) {
        self.viewer.hide(structure.name);
        self.viewer.requestRedraw();
        
        if(structure == self.reference && self.structures.length > 0) {
            self.reference = null;
        }
        
        return structure;
    }
    
    self.removeStructure = function(structure) {
        self.viewer.rm(structure.name);
        self.viewer.requestRedraw()
        
        self.structures = self.structures.filter(function (el) {
            return el.name !== structure.name;
        });
        
        structure.struc = null;
        structure.geom = null;
        
        return structure;
    }
    
    self.getStructureByAtom = function(atom) {
        for(var i = 0; i < self.structures.length; i++) {
            var s = self.structures[i];
            if(s.geom != null) {
            
                if(s.geom.structure() == atom.structure()) { 
                    return  s;
                }
            }
        }
        return null;
    }
    
    self.loadViewer = function(refresh) {
        if(self.viewer == null || refresh == true) {
            var options = {
                width: self.options.width,
                height: self.options.struc_height,
                antialias: true,
                quality : self.options.struc_quality,
                animateTime: 1000
            };
            
            var parent = self.viewer_container;
            parent.innerHTML = '';
            
            self.viewer = pv.Viewer(parent, options);
            
            // highlight residues under the cursor
            parent.addEventListener('mousemove', function(event) {
                //get the atom
                var rect = self.viewer.boundingClientRect();
                var picked = self.viewer.pick({ x : event.clientX - rect.left, y : event.clientY - rect.top });
                
                //if it is the same atom as previously, don't do anything
                if (self.prevPicked !== null && picked !== null &&
                    picked.target() === self.prevPicked.atom) {
                    return;
                }
                
                // reset color of previously picked atom.
                if (self.prevPicked !== null) {
                  self.setColorForAtom(self.prevPicked.node, self.prevPicked.atom, self.prevPicked.color);
                }
                
                if (picked !== null) {
                    var atom = picked.target();
                    var s = self.getStructureByAtom(atom);
                  
                    var atom_name = atom.qualifiedName();
                    console.log(atom_name)
                    console.log(atom_name[0])
                    if(atom_name[0] == " ")
                        atom_name = atom_name.substring(2);
                  
                    self.atom_name.innerHTML = s.name + ' - ' + atom_name;
                  
                    // store original color before changing it to the highlight color.
                    var color = [0,0,0,0];
                    picked.node().getColorForAtom(atom, color);
                    self.prevPicked = { atom : atom, color : color, node : picked.node() };
            
                    self.setColorForAtom(picked.node(), atom, 'red');
                } else {
                    self.atom_name.innerHTML = '&nbsp;';
                    self.prevPicked = null;
                }
                
                self.viewer.requestRedraw();
            });
            
            // select clicked residues
            self.viewer.on('click', function(picked, ev) {
                //do nothing if no atom was selected
                if (picked === null || picked.target() === null) {
                    return;
                }
                
                var atom = picked.target();
                var struc = self.getStructureByAtom(atom);
                
                if (struc == null) {
                    return;
                }
                
                var sel = struc.geom.selection();
                
                // get residue position
                var struc_pos = atom.residue().index();
                var seq_pos = self.getAlignPosFromStrucPos(struc, struc_pos);
                
                if (!sel.removeAtom(atom, true)) {
                    self.selectResidue(seq_pos);
                } else {
                    self.deselectResidue(seq_pos);
                }
            });
        }
    }
    
    self.setColorForAtom = function(go, atom, color) {
        var view = go.structure().createEmptyView();
        view.addAtom(atom);
        go.colorBy(pv.color.uniform(color), view);
    }
    
    self.loadMSA = function() {
        self.alignment = new msa({ el: self.alignment_container });
        
        self.alignment.g.zoomer.set("alignmentHeight", self.options.seq_height); 
        self.alignment.g.zoomer.set("labelIdLength", 0); 
        self.alignment.g.zoomer.set("labelNameLength", 100); 
        self.alignment.g.zoomer.set("residueFont", 15); 
        self.alignment.g.zoomer.set("labelFontsize", 14); 
        self.alignment.g.zoomer.set("canvasEventScale", 2); 
        
        self.alignment.g.vis.set("labelId", false);
        
        self.alignment.g.config.set("registerMouseHover", true);
        self.alignment.g.on("column:click", self.selectAlignmentResidues);
        self.alignment.g.on("row:click", self.selectAlignmentResidues);
        self.alignment.g.on("residue:click", self.selectAlignmentResidues);
        
        self.alignment.render();
    }
    
    self.selectAlignmentResidues = function(col) {
        var pos = col.rowPos
        if(self.columns[pos] == true) {
            self.deselectResidue(pos);
        } else {
            self.selectResidue(pos)
        }
        
        self.updateAlignmentSelection();
    }
    
    self.getAlignPosFromStrucPos = function(structure, pos) {
        for(var i = pos; i < structure.residues.length; i++) {
            var residue = structure.residues[i];
            
            if(residue.pos == pos + 1) {
                return i;
            }
        }
    }
    
    self.getStrucPosFromAlignPos = function(structure, pos) {
        var residues = structure.residues;
        return residues[pos].pos;
    }
    
    self.selectResidue = function(index) {
        self.columns[index] = true;
        
        //select structure residues
        for(var j = 0; j < self.structures.length; j++) {
            var struc = self.structures[j];
            if(struc.geom != null) {
                self.setResidueSelected(struc, index, true);
            }
        }
        
        self.viewer.requestRedraw();
        
        self.updateAlignmentSelection();
    }
    
    self.deselectResidue = function(index) {
        self.columns[index] = false;
        
        //select structure residues
        for(var j = 0; j < self.structures.length; j++) {
            var struc = self.structures[j];
            if(struc.geom != null) {
                self.setResidueSelected(struc, index, false);
            }
        }
        
        self.viewer.requestRedraw();
        
        self.updateAlignmentSelection();
    }
    
    self.setResidueSelected = function(struc, index, selected) {
        var sel = struc.geom.selection();
                
        var pos = self.getStrucPosFromAlignPos(struc, index);
        if(pos > -1) {
            if(selected) {
                var subset = struc.struc.select({ rindexRange: [pos - 1, pos - 1] });
            } else {
                var subset = sel.select({ rnumRange: [pos, pos] });
            }
                
            for(var i = 0; i < subset.atoms().length; i++) { 
                var atom = subset.atoms()[i];
                
                if(selected) {
                    sel.addAtom(atom);
                } else {
                    sel.removeAtom(atom, true);
                }
            }
            
            struc.geom.setSelection(sel);
        }
    }
    
    self.updateAlignmentSelection = function() {
        self.alignment.g.selcol.reset();
        for(var i = 0; i < self.columns.length; i++) {
            if(self.columns[i] == true) {
                self.alignment.g.selcol.add(
                    new msa.selection.columnsel({ xStart: i, xEnd: i })
                );
            }
        }
    }
    
    self.showLoading = function() {
        self.viewer_container.className = "masked"
    }
    
    self.hideLoading = function() {
        self.viewer_container.className = ""
    }
}
