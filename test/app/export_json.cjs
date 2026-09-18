const assert = require("node:assert/strict");
const fs = require("node:fs");
const vm = require("node:vm");
const test = require("node:test");

test("export contains integer indices, exact lattice geometry, and edges", () => {
  const context = vm.createContext({ document: { querySelector: () => ({}) } });
  const source = fs.readFileSync(`${__dirname}/../../app/app.js`, "utf8");
  // Load the editor functions without registering browser events or rendering.
  vm.runInContext(source.slice(0, source.indexOf('\ndocument.querySelectorAll(".tool")')), context);
  vm.runInContext(`
    state.nodes = new Map([
      ["0:0", { id: "0:0", q: 0, r: 0, weight: 1 }],
      ["0:1", { id: "0:1", q: 0, r: 1, weight: 2, x: 99, y: 99 }],
      ["1:1", { id: "1:1", q: 1, r: 1, weight: 1 }],
      ["0:2", { id: "0:2", q: 0, r: 2, weight: 1 }],
    ]);
    state.pins = ["0:2", "0:0"];
  `, context);
  const exported = () => JSON.parse(vm.runInContext("JSON.stringify(payload())", context));
  const triangular = exported();
  assert.deepEqual(triangular.nodes, [
    { id: "0:0", q: 0, r: 0, weight: 1 },
    { id: "0:1", q: 0, r: 1, weight: 2 },
    { id: "1:1", q: 1, r: 1, weight: 1 },
    { id: "0:2", q: 0, r: 2, weight: 1 },
  ]);
  assert.equal(triangular.lattice.index_base, 0);
  assert.deepEqual(triangular.lattice.basis, [[1, 0], [0, "sqrt(3)/2"]]);
  assert.deepEqual(triangular.lattice.odd_row_offset, ["1/2", 0]);
  assert.deepEqual(triangular.edges, [
    { source: "0:0", target: "0:1" },
    { source: "0:1", target: "1:1" },
    { source: "0:1", target: "0:2" },
  ]);
  assert.deepEqual(triangular.pins, ["0:2", "0:0"]);
  vm.runInContext('state.latticeShape = "KSG"', context);
  const square = exported();
  assert.deepEqual(square.nodes, triangular.nodes);
  assert.deepEqual(square.lattice.basis, [[1, 0], [0, 1]]);
  assert.deepEqual(square.lattice.odd_row_offset, [0, 0]);
  assert.equal(square.edges.length, 5);
});
