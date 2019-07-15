import * as systre from './src/scripts/systreCmd';

import { Archive } from './src/io/archive';
import * as symmetries from './src/pgraphs/symmetries';
import * as periodic from './src/pgraphs/periodic';
import parseDSymbols from './src/io/ds';
import * as tilings from './src/dsymbols/tilings';
import * as cgd from './src/io/cgd';
import { systreKey } from './src/pgraphs/invariant';

const NEW_TOPOLOGY_CODE = 'UNKNOWN';  // formerly 'NEW'
const MISMATCH_TOPOLOGY_CODE = 'MISMATCH';

// One-liner from https://stackoverflow.com/questions/14832603/check-if-all-values-of-array-are-equal
const allEqual = arr => arr.every( v => v === arr[0] );

const processRCSR = (input, options, archives = []) => {
  // Re-implementing systre.processGraph by extracting the relevant data
  const { graph, name, nodes: originalNodes } = input;
  if (!checkGraph(graph, console.log)) {
    return "ERROR";
  }
  const { graph: G, orbits: translationOrbits } =
    symmetries.minimalImageWithOrbits(graph);

  const key = systreKey(G);
  //const countMatches = showAndCountGraphMatches(key, archives, writeInfo);
  const found = archives[0].getByKey(key);
  if (found) {
    return found.id;
  } else {  // 'Structure is new for this run.'
    // archives.find(arc => arc.name == '__internal__').addNet(G, name, key);  // from systre.processGraph
    archives[0].addNet(G, NEW_TOPOLOGY_CODE, key);
    return NEW_TOPOLOGY_CODE;
  }
};


export const getRCSRTopology = (topology_cgd, rcsr_archive) => {
  // Extract the RCSR topology, similar to systre.processData
  const archive_obj =  new Archive('RCSRnets.arc');
  archive_obj.addAll(rcsr_archive);
  const archives = [archive_obj, new Archive('__internal__')];

  const options = {
    outputEmbedding: true,
    relaxPositions: true,
    skipWarnings: true
  };

  // systre.processData(topology_cgd, "topology.cgd", options, archives);
  const inputs = nets(topology_cgd, "topology.cgd");
  let topologies = [];
  for (const input of inputs) {
    if (input.errors.length != 0) {
      return "ERROR";
    }
    try {
      //const process = periodic.isConnected(input.graph) ? processGraph : processDisconnectedGraph;
      let components = [];
      if (periodic.isConnected(input.graph)) {
        components = [{graph: input.graph}];
      } else {  // from systre.processDisconnectedGraph
        components = periodic.connectedComponents(input.graph);
      }
      for (let i = 1; i <= components.length; ++i) {
        const comp = components[i - 1];  // Using 1-indexing instead of 0-indexing for names below
        topologies.push(processRCSR(
          {graph: comp.graph, name: `${input.graph.name}_component_${i}`},
          options, archives
        ));
      }
    } catch(ex) {
      systre.reportSystreError('INTERNAL', ex + '\n' + ex.stack, console.log);
      return "ERROR";
    }

    if (allEqual(topologies)) {
      return topologies[0];
    } else {
      return MISMATCH_TOPOLOGY_CODE;
    }
  }
};

const nets = function*(data, fileName) {
  // Copied from systre.nets
  if (fileName.match(/\.(ds|tgs)$/))
    for (const t of parseDSymbols(data)) {
      const g = tilings.skeleton(tilings.makeCover(t.symbol));
      yield Object.assign({ warnings: [], errors: [] }, g, t);
    }
  else if (fileName.match(/\.(cgd|pgr)$/))
    for (const g of cgd.structures(data)) {
      yield g;
    }
};

const checkGraph = (graph, writeInfo) => {
  // Copied from systre.checkGraph
  if (!periodic.isLocallyStable(graph)) {
    const msg = ("Structure has collisions between next-nearest neighbors."
                 + " Systre does not currently support such structures.");
    reportSystreError("STRUCTURE", msg, writeInfo);
    return false;
  }

  if (symmetries.isLadder(graph)) {
    const msg = "Structure is non-crystallographic (a 'ladder')";
    reportSystreError("STRUCTURE", msg, writeInfo);
    return false;
  }

  if (periodic.hasSecondOrderCollisions(graph)) {
    const msg = ("Structure has second-order collisions."
                 + " Systre does not currently support such structures.");
    reportSystreError("STRUCTURE", msg, writeInfo);
    return false;
  }

  if (!periodic.isStable(graph)) {
    writeInfo("Structure has collisions.");
    writeInfo();
  }

  return true;
};

const reportSystreError = (errorType, message, writeInfo) => {
  // Copied from systre.reportSystreError
  writeInfo("==================================================");
  writeInfo(`!!! ERROR (${errorType}) - ${message}`);
  writeInfo("==================================================");
  writeInfo();
};

