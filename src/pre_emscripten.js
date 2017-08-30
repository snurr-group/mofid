/* Set up environment variables */
/* https://kripken.github.io/emscripten-site/docs/porting/connecting_cpp_and_javascript/Interacting-with-code.html#interacting-with-code-environment-variables */
Module.preRun.push(function() {ENV.BABEL_DATADIR = "/ob_datadir"})
