import * as systre from './src/scripts/systreCmd';

const getRCSRTopology = (topology_cgd, rcsr_archive) => {
  let archives = [rcsr_archive];

  const options = {
    outputEmbedding: true,
    relaxPositions: true,
    skipWarnings: true
  };

  systre.processData(topology_cgd, "topology.cgd", options, archives);
  // TODO: save the info to a variable with two new methods writing strings to global variables instead of
  // using the default prefixedLineWriter to console.log

  return "STUB";
};

// Export the function to the global namespace
// This feature could also be implemented via the "library" function, but this is easier
// See also https://stackoverflow.com/questions/33148787/how-expose-a-exported-function-into-global-scope-with-babel-and-webpack
window.getRCSRTopology = getRCSRTopology;
// WARNING: it might not work properly due to exporting the proper functions.
// Let's check for sure using a build-dev to figure out where things crash
// And maybe look at the proper implementation method: https://webpack.js.org/configuration/output/#expose-a-variable

