import * as systre from './src/scripts/systreCmd';
import { Archive } from './src/io/archive';

export const getRCSRTopology = (topology_cgd, rcsr_archive) => {
  const archive_obj =  new Archive('RCSRnets.arc');
  archive_obj.addAll(rcsr_archive);
  const archives = [archive_obj, new Archive('__internal__')];

  const options = {
    outputEmbedding: true,
    relaxPositions: true,
    skipWarnings: true
  };

  systre.processData(topology_cgd, "topology.cgd", options, archives);
  // TODO: save the info to a variable with two new methods writing strings to global variables instead of
  // using the default prefixedLineWriter to console.log - nah, just define the overwriting function here
  // And probably two new global functions of getRCSRTopology and runSystre
  // Also, need to actually use the topology (other than just printing the Systre results -- start by saving a string and printing it to the thing)

  return "STUB";
};

