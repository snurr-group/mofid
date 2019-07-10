var fs = require('fs');
var path = require('path');

const webpack = require('webpack');
const CleanWebpackPlugin = require('clean-webpack-plugin');

var execSync = require('child_process').execSync;
var execOpts = { encoding: 'utf8' };

// Remove git dependency (because we don't want nested .git repos)
// by including the metadata from /Resources/webGavrog-20190709.zip
//var gitRev = execSync('git rev-parse HEAD', execOpts);
var gitRev = "f3cfb3082de510b1b7f319f9ab28429b08c038db";
//var gitDate = execSync('git show -s --format=%ci HEAD', execOpts);
var gitDate = "2019-07-02 17:58:15 +1000";


fs.writeFileSync(path.resolve(__dirname, 'src', 'version.js'),
                 'export const gitRev = "' + gitRev.trim() + '";\n' +
                 'export const gitDate = "' + gitDate.trim() + '";\n');


module.exports = {
  entry: ['babel-polyfill', './main.js'],
  node: { fs: 'empty' },  // per https://github.com/webpack-contrib/css-loader/issues/447 via Google
  plugins: [
    new CleanWebpackPlugin()
  ],
  output: {  // export as a library, per https://webpack.js.org/guides/author-libraries/
    filename: '[name].js',
    path: path.resolve(__dirname, 'dist'),
    library: 'webGavrog',
    libraryTarget: 'var'
  },
  optimization: {
    runtimeChunk: 'single',
    splitChunks: {
      cacheGroups: {
        vendor: {
          test: /[\\/]node_modules[\\/]/,
          name: 'vendors',
          chunks: 'all'
        }
      }
    }
  },
  module: {
    rules: [
      { test: /\.js$/,
        exclude: /node_modules/,
        use: "babel-loader"
      },
      { test: /\.elm$/,
        exclude: [/node_modules/, /elm-stuff/],
        use: "elm-webpack-loader"
      },
      {
        test: /sceneWorker\.js$/,
        use: {
          loader: 'worker-loader',
          options: {
            inline: true
          }
        }
      }
    ]
  },
  resolve: {
    extensions: [ ".js", ".elm" ]
  }
}
