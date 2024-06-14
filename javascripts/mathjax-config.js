window.MathJax = {
    tex: {
        tags: 'ams',
        inlineMath: [['$', '$'], ['\\(', '\\)']],
        displayMath: [['$$', '$$'], ['\\[', '\\]']],
        macros: {
            Re: "\\operatorname{Re}",
            Im: "\\operatorname{Im}"
        }
    }
};
  
/* -------------------------------------------------------------------------- */

// window.MathJax = {
//     tex: {
//         tags: 'ams',
//         inlineMath: [['$', '$'], ['\\(', '\\)']],
//         displayMath: [['$$', '$$'], ['\\[', '\\]']],
//         // processEscapes: true,
//         // processEnvironments: true,
//         // processRefs: true,
//         // packages: {'[+]': ['noerrors']}
//     },
//     options: {
//         // ignoreHtmlClass: 'tex2jax_ignore',
//         // processHtmlClass: 'tex2jax_process',
//         // ignoreHtmlClass: ".*|",
//         // processHtmlClass: "arithmatex"
//     },
//     // loader: {load: ['[tex]/noerrors']}
// };


// // document$.subscribe(() => { 
// //     MathJax.startup.output.clearCache()
// //     MathJax.typesetClear()
// //     MathJax.texReset()
// //     MathJax.typesetPromise()
// //   })


// // document$.subscribe(({ body }) => { 
// // renderMathInElement(body, {
// //     delimiters: [
// //     { left: "$$",  right: "$$",  display: true },
// //     { left: "$",   right: "$",   display: false },
// //     { left: "\\(", right: "\\)", display: false },
// //     { left: "\\[", right: "\\]", display: true }
// //     ],
// // })
// // })