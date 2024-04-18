<script>
  import { onMount, afterUpdate } from "svelte";
  import * as d3 from "d3";
  import * as math from "mathjs";
  import Toggle from "./Toggle.svelte";
  import ft from "fourier-transform";

  import {
    linspace,
    Sellimeier_PPLN,
    Sellimeier_PPKTP,
    AiryFunction,
    finesse,
    WL_by_energy_conservation,
    create_50ghz_wl,
    n_raicol_ppktp,
    sellimeier_ppln_single,
    n_ppln,
    phase_matching,
    RaicolNDataPPKTP,
    FilterComb,
    interp
  } from "./util.js";

  let canvas; // reference to the canvas element
  let svg; // reference to the svg element
  let cvs_width;
  let cvs_height;
  let isMounted = false;

  let signalCDS;
  let idlerCDS;
  let dwdm_filterCDS;
  let indexCDS;
  let fourierCDS;

  let backup_array;
  let old_show_2d = true;
  let previousParams = {};
  let canvas_initialized = false;
  let dwdm_filter_spectrum;
  let n_raicol_ppktp_data;
  let index_function;

  let filter_comb_idler;
  let filter_comb_signal;

  let raicol_n;

  const params = {
    show_2d: true,
    c: 299792458, // speed of light in m/s
    WLpumpair: 774.08,
    numWLs: 501,
    BW: 10e-9, // bandwidth in m
    WLSigma: 0.815, // sigma in nm
    L: 1, // length in mm
    L_multiplier: 1.62, // 1.71,
    R1: 0.85,
    R2: 0.1,
    Texpt: 88, //24,
    pump_checked: true,
    phase_matching_checked: true,
    waveguide_checked: true,
    ppktp_checked: true,
    start_idler: 1548.51,
    end_idler: 1568.36,
    start_signal: 1528.77,
    end_signal: 1548.11,
    polling_period: 24.838
  };

  const wl_50ghz = create_50ghz_wl(); //goes from short to long wl in nm
  // console.log("wl_50ghz: ", wl_50ghz)

  const WLCenteral = params.WLpumpair * 1e-9 * 2;

  function compute_result(p, init) {
    const Omega0 = (2 * Math.PI * p.c) / (2 * p.WLpumpair * 1e-9);
    const WLCenteral = p.WLpumpair * 1e-9 * 2;

    if (isNaN(p.L_multiplier)) {
      return 1;
    }

    const polling_period = p.polling_period*1e-6;

    const WLSigma = p.WLSigma * 1e-9;
    const L = p.L * 1e-3 * p.L_multiplier;

    const WLSignalairJSI = linspace(
      p.start_signal * 1e-9,
      p.end_signal * 1e-9,
      p.numWLs,
    );

    const WLIdlerairJSI = linspace(
      p.start_idler * 1e-9,
      p.end_idler * 1e-9,
      p.numWLs,
    );

    let nSignal = new Array(WLSignalairJSI.length).fill(1);
    let nIdler = new Array(WLSignalairJSI.length).fill(1);

    if (p.ppktp_checked) {
      index_function = raicol_n.index_function;
    } else {
      index_function = n_ppln;
    }

    nSignal = WLSignalairJSI.map((wl) => index_function(wl, p.Texpt));
    nIdler = WLIdlerairJSI.map((wl) => index_function(wl, p.Texpt));

    const ASignal = AiryFunction(p.R1, p.R2, L, nSignal, WLSignalairJSI);
    const AIdler = AiryFunction(p.R1, p.R2, L, nIdler, WLIdlerairJSI);

    const OmegaSignalJSI = WLSignalairJSI.map((wl) => (2 * Math.PI * p.c) / wl);
    let OmegaIdlerJSI = WLIdlerairJSI.map((wl) => (2 * Math.PI * p.c) / wl);
    const OmegaSignalGrid = math.matrix(OmegaSignalJSI);
    const OmegaIdlerGrid = math.matrix(OmegaIdlerJSI);

    // const p.
    const FSigma = (WLSigma * p.c) / (2 * p.WLpumpair * 1e-9) ** 2;

    // const full_spectrum = linspace(
    //   wl_50ghz[0],
    //   wl_50ghz[wl_50ghz.length - 1],
    //   501,
    // );

    const full_spectrum = linspace(500, 3500, 501);
    const full_index = full_spectrum.map((value) =>
      index_function(value * 1e-9, p.Texpt, false),
    );

    // console.log("full index: ", full_index)
    // console.log(full_spectrum.map((value) => value * 1e-9))
    // console.log("index: ", full_index)
    // console.log("full spectrum: ", full_spectrum)
    // console.log("index T: ", p.Texpt)

    if (
      signalCDS?.data?.x &&
      !init &&
      signalCDS.data.x.length == p.numWLs &&
      old_show_2d == p.show_2d
    ) {
      const filter_comb_idler_and_airy = filter_comb_idler.out_array.map(
        (value, i) => value * AIdler[i],
      );
      const filter_comb_signal_and_airy = filter_comb_signal.out_array.map(
        (value, i) => value * ASignal[i],
      );

      indexCDS.data.x = full_spectrum;
      indexCDS.data.y = full_index;
      signalCDS.data.y = ASignal;
      idlerCDS.data.y = AIdler;

      dwdm_filterCDS.data.y_idler_airy = filter_comb_idler_and_airy;
      dwdm_filterCDS.data.y_signal_airy = filter_comb_signal_and_airy;
      dwdm_filterCDS.data.asignal = ASignal;
      dwdm_filterCDS.data.aidler = AIdler;

      const wlSection = WLIdlerairJSI.slice(Math.floor(p.numWLs*.01),Math.floor(p.numWLs*.0306))
      const transSection = filter_comb_idler_and_airy.slice(Math.floor(p.numWLs*.01),Math.floor(p.numWLs*.0306))

      // const transF = ft(transSection)

      const x_256 = linspace(wlSection[0], wlSection[wlSection.length - 1], 256);
      const y_256 = interp(x_256, wlSection, transSection);


      const transF = ft(y_256)

      const x_256_f = x_256.map((value) => p.c/value);
      const df = x_256_f[x_256_f.length - 2] - x_256_f[x_256_f.length - 1]
      // console.log("df: ", df)
      const dt = 1/(df*255)
      const t = linspace(0, dt*256, 256).map((value) => value * 1e12)

      fourierCDS.data.x = t.slice(0,128);
      // fourierCDS.data.y = transSection;
      fourierCDS.data.y = transF;

      dwdm_filterCDS.change.emit();
      signalCDS.change.emit();
      idlerCDS.change.emit();
      indexCDS.change.emit();
      fourierCDS.change.emit();
    } else {
      filter_comb_idler = new FilterComb(
        dwdm_filter_spectrum,
        WLIdlerairJSI.map((value) => value * 1e9),
      );
      filter_comb_signal = new FilterComb(
        dwdm_filter_spectrum,
        WLSignalairJSI.map((value) => value * 1e9),
      );

      const filter_comb_idler_and_airy = filter_comb_idler.out_array.map(
        (value, i) => value * AIdler[i],
      );

      const filter_comb_signal_and_airy = filter_comb_signal.out_array.map(
        (value, i) => value * ASignal[i],
      );

      const wlSection = WLIdlerairJSI.slice(Math.floor(p.numWLs*.01),Math.floor(p.numWLs*.0306))
      const transSection = filter_comb_idler_and_airy.slice(Math.floor(p.numWLs*.01),Math.floor(p.numWLs*.0306))

      const x_256 = linspace(wlSection[0], wlSection[wlSection.length - 1], 256);
      const y_256 = interp(x_256, wlSection, transSection);


      const transF = ft(y_256)

      const x_256_f = x_256.map((value) => p.c/value);
      const df = x_256_f[x_256_f.length - 2] - x_256_f[x_256_f.length - 1]
      // console.log("df: ", df)
      const dt = 1/(df*255)
      const t = linspace(0, dt*256, 256).map((value) => value * 1e12)



      if (!p.show_2d) {
        const inner_off_array = Array(OmegaIdlerJSI.length).fill(1);
        backup_array = Array(OmegaSignalJSI.length).fill(inner_off_array);
      }
      old_show_2d = p.show_2d;

      signalCDS = new Bokeh.ColumnDataSource({
        data: {
          x: WLSignalairJSI.map((value) => value * 1e9),
          y: ASignal,
        },
      });

      idlerCDS = new Bokeh.ColumnDataSource({
        data: {
          x: WLIdlerairJSI.map((value) => value * 1e9),
          y: AIdler,
        },
      });

      // dwdm_filterCDS = new Bokeh.ColumnDataSource({
      //   data: {
      //     x: dwdm_filter_spectrum.x,
      //     y: dwdm_filter_spectrum.y,
      //   },
      // });
      dwdm_filterCDS = new Bokeh.ColumnDataSource({
        data: {
          x_idler: WLIdlerairJSI.map((value) => value * 1e9),
          x_signal: WLSignalairJSI.map((value) => value * 1e9),
          y_idler: filter_comb_idler.out_array,
          y_signal: filter_comb_signal.out_array,

          y_idler_airy: filter_comb_idler_and_airy,
          y_signal_airy: filter_comb_signal_and_airy,

          asignal: ASignal,
          aidler: AIdler,
        },
      });

      indexCDS = new Bokeh.ColumnDataSource({
        data: {
          x: full_spectrum,
          y: full_index,
        },
      });

      fourierCDS = new Bokeh.ColumnDataSource({
        data: {
          // x: wlSection,
          // y: transSection,
          x: t.slice(0,128),
          y: transF,
          // y: transF,
        },
      });

      // const x_range = [wl_50ghz[0], wl_50ghz[wl_50ghz.length - 1]];

      const x_range = new Bokeh.Range1d({
        start: wl_50ghz[0],
        end: wl_50ghz[wl_50ghz.length - 1],
      });

      // const range_thing = new Bokeh.Plotting.Models.range1d({
      //   start: wl_50ghz[0],
      //   end: wl_50ghz[wl_50ghz.length - 1],
      // });
      const scope = new Bokeh.Plotting.figure({
        // title: "Example of random data",
        // tools: "pan,wheel_zoom,box_zoom,reset,save",
        tools: "xwheel_zoom, xpan, reset",
        // toolbar: {logo: null},
        //tools: "",
        // sizing_mode: "stretch_both",
        active_drag: "xwheel_zoom",
        sizing_mode: "stretch_width",
        height: 250,
        x_range: x_range,
        y_range: [0, d3.max(ASignal)],
        output_backend: "webgl",
      });

      const filter_graph = new Bokeh.Plotting.figure({
        tools: "xwheel_zoom, xpan, reset",
        active_drag: "xwheel_zoom",
        sizing_mode: "stretch_width",
        height: 250,
        x_range: x_range,
        y_range: [0, 1],
        output_backend: "webgl",
      });

      // const index_graph = new Bokeh.Plotting.figure({
      //   height: 300,
      //   x_range: [wl_50ghz[0], wl_50ghz[wl_50ghz.length - 1]],
      //   y_range: [0, 5],
      // });

      const index_graph = new Bokeh.Plotting.figure({
        title: "Index of Refraction nz",
        height: 300,
        x_range: [500, 3500],
        y_range: [1.75, 1.925],
        x_axis_label: "Wavelength (nm)",
        y_axis_label: "Index of Refraction",
      });

      const fourier_graph = new Bokeh.Plotting.figure({
        title: "Half of pulse in time-domain",
        height: 300,
        x_range: [0, 500],
        // y_range: [1.65, 1.825],
        x_axis_label: "Time (ps)",
        y_axis_label: "Intensity",
      });

      scope.toolbar.logo = null;
      const line_1 = scope.line(
        { field: "x" },
        { field: "y" },
        {
          source: signalCDS,
          line_width: 3,
          line_color: "#333333",
        },
      );

      const line_2 = scope.line(
        { field: "x" },
        { field: "y" },
        {
          source: idlerCDS,
          line_width: 3,
          line_color: "#330033",
        },
      );

      const filter_line_1 = filter_graph.line(
        { field: "x_idler" },
        { field: "y_idler" },
        {
          line_width: 3,
          line_alpha: 0.3,
          source: dwdm_filterCDS,
          line_color: "#FF3333",
        },
      );

      const filter_line_2 = filter_graph.line(
        { field: "x_signal" },
        { field: "y_signal" },
        {
          line_width: 3,
          line_alpha: 0.3,
          source: dwdm_filterCDS,
          line_color: "#FF3333",
        },
      );

      const filter_line_passage_1 = filter_graph.line(
        { field: "x_signal" },
        { field: "y_signal_airy" },
        {
          line_width: 3,
          source: dwdm_filterCDS,
          line_color: "#3333FF",
        },
      );

      const filter_line_passage_2 = filter_graph.line(
        { field: "x_idler" },
        { field: "y_idler_airy" },
        {
          line_width: 3,
          source: dwdm_filterCDS,
          line_color: "#33FF33",
        },
      );

      const filter_line_asingal = filter_graph.line(
        { field: "x_signal" },
        { field: "asignal" },
        {
          line_width: 3,
          line_alpha: 0.3,
          source: dwdm_filterCDS,
          line_color: "#333333",
        },
      );

      const filter_line_aidler = filter_graph.line(
        { field: "x_idler" },
        { field: "aidler" },
        {
          line_width: 3,
          line_alpha: 0.3,
          source: dwdm_filterCDS,
          line_color: "#333333",
        },
      );

      const index_line = index_graph.line(
        { field: "x" },
        { field: "y" },
        {
          source: indexCDS,
          line_width: 3,
          line_color: "#333333",
        },
      );

      const trans_line = fourier_graph.line(
        { field: "x" },
        { field: "y" },
        {
          source: fourierCDS,
          line_width: 3,
          line_color: "#333333",
        },
      );

      const doc = new Bokeh.Document();
      doc.add_root(scope);
      doc.add_root(filter_graph);
      doc.add_root(index_graph);
      doc.add_root(fourier_graph);

      const plotElement = document.getElementById("plot");

      // Check if plotElement exists before trying to remove its child nodes
      if (plotElement) {
        while (plotElement.firstChild) {
          plotElement.firstChild.remove();
        }
      }

      Bokeh.embed.add_document_standalone(doc, plotElement);
    }

    // Calculate the JSI using the Gaussian approximation

    // console.log("wlSignalairJSI: ", WLSignalairJSI)
    // WLSignalairJSI is in meters

    // this is important!!!
    // the frequency array got flipped, but nothing else did!!!
    OmegaIdlerJSI = OmegaIdlerJSI.reverse();
    let JSI = p.show_2d
      ? OmegaSignalJSI.map((os, i) => {
          return OmegaIdlerJSI.map((oi, j) => {
            let debug_show = false;
            if (i == 100 && j == 501 - 100) {
              debug_show = true;
            }

            if (i == 400 && j == 501 - 400) {
              debug_show = true;
            }
            // if (Math.abs(i - 100) <= 5 && Math.abs(j - (501 - 100)) <= 5) {
            //   debug_show = true;
            // }

            // if (Math.abs(i - 400) <= 5 && Math.abs(j - (501 - 400)) <= 5) {
            //   debug_show = true;
            // }

            const waveguide = p.waveguide_checked ? AIdler[j] * ASignal[i] : 1;
            const pump = p.pump_checked
              ? Math.exp(
                  -Math.pow(os + oi - 2 * Omega0, 2) /
                    (2 * Math.pow(FSigma, 2)),
                )
              : 1;
            const ps = p.phase_matching_checked
              ? phase_matching(
                  (2 * Math.PI * p.c) / oi,
                  (2 * Math.PI * p.c) / os,
                  index_function,
                  p.Texpt,
                  polling_period,
                  L,
                  debug_show,
                )
              : 1;
            // if (i == 200 && j == 200) {
            //   console.log("waveguide: ", waveguide)
            //   console.log("pump: ", pump)
            //   console.log("ps: ", ps)
            //   console.log("index function: ", index_function(1550e-9, p.Texpt))
            //   console.log("n_ppln: ", sellimeier_ppln_single(1550e-9, "eray", p.Texpt))
            //   console.log("ppktp: ", n_raicol_ppktp(1550e-9, p.Texpt))
            //   console.log("p.Texpt: ", p.Texpt)
            //   console.log("sellimeier_ppln_single", sellimeier_ppln_single(1550e-9, "eray", p.Texpt))
            // }
            return waveguide * pump * ps;
          });
        })
      : backup_array;

    if (p.show_2d) {
      // still trying to figure out how to initialize 2d plot
      

      // see here for syntax that last worked

      // const source_im = new Bokeh.ColumnDataSource({
      //   data: {
      //     image: JSI,
      //   },
      // });

      // console.log(Bokeh);

      

      // const source = new Bokeh.ColumnDataSource({
      //   data: {
      //     image: [
      //       [
      //         [0, 1, 0.25],
      //         [1, 0, 0.75],
      //       ],
      //     ],
      //     x: [0],
      //     x2: [1],
      //     x3: [2],
      //     y: [1],
      //     dw: [0.8],
      //     dh: [1],
      //   },
      // });

      // var plott = Bokeh.Plotting.figure({
      //   title: "Example of Random data",
      //   tools: "pan,wheel_zoom,box_zoom,reset,save",
      //   height: 300,
      //   width: 900,
      // });

      // // add a line with data from the source
      // plott.image({
      //   image: { field: "image" },
      //   x: { field: "x" },
      //   y: { field: "y" },
      //   dw: { field: "dw" },
      //   dh: { field: "dh" },
      //   // color_mapper: cmap,
      //   source: source,
      // });

      // const doc = new Bokeh.Document();
      // doc.add_root(plott);
      // const plotElement = document.getElementById("plot_2");
      // Bokeh.embed.add_document_standalone(doc, plotElement);

      canvas.style.opacity = 1;
      imshow(
        JSI,
        1 * (501 / p.numWLs),
        d3.scaleSequential(d3.interpolateViridis),
      );
    } else {
      canvas.style.opacity = 0.2;
    }
  }

  $: if (isMounted && dwdm_filter_spectrum) {
    if (JSON.stringify(params) !== JSON.stringify(previousParams)) {
      // console.log('params changed');
      previousParams = { ...params };
      compute_result(params, false);
    }
  }

  onMount(async () => {
    const filter_response = await fetch("./50_GHz_spectrum.json");
    dwdm_filter_spectrum = await filter_response.json();

    const n_ppktp_data = await fetch("./ppKTP_n.json").then((response) =>
      response.json(),
    );

    raicol_n = new RaicolNDataPPKTP(n_ppktp_data);

    compute_result(params, true);

    // Create axes
    const xAxis = d3.axisBottom(
      d3
        .scaleLinear()
        .domain([params.start_signal, params.end_signal])
        .range([0, cvs_width]),
    );
    const yAxis = d3.axisLeft(
      d3
        .scaleLinear()
        .domain([params.start_idler, params.end_idler])
        .range([cvs_height, 0]),
    );

    // // Append axes to SVG
    // d3.select(svg)
    //   .append("g")
    //   .attr("transform", `translate(100,${100 + cvs_height})`)
    //   .call(xAxis);

    // d3.select(svg)
    //   .append("g")
    //   .attr("transform", `translate(100,100)`)
    //   .call(yAxis);
    // Append axes to SVG
    const xAxisGroup = d3
      .select(svg)
      .append("g")
      .attr("transform", `translate(100,${100 + cvs_height})`)
      .call(xAxis);

    const yAxisGroup = d3
      .select(svg)
      .append("g")
      .attr("transform", `translate(100,100)`)
      .call(yAxis);

    // Append labels to the axes
    xAxisGroup
      .append("text")
      .attr("class", "axis-label")
      .attr("x", cvs_width / 2)
      .attr("y", 40)
      .style("text-anchor", "middle")
      .style("fill", "black") // Ensure the text is black
      .text("Signal (nm)");

    yAxisGroup
      .append("text")
      .attr("class", "axis-label")
      .attr("x", -cvs_height / 2)
      .attr("y", -60)
      .attr("transform", "rotate(-90)")
      .style("text-anchor", "middle")
      .style("fill", "black") // Ensure the text is black
      .text("Idler (nm)");

    // Set the width and height of the SVG
    svg.setAttribute("width", `${cvs_width + 200}px`);
    svg.setAttribute("height", `${cvs_height + 200}px`);
    isMounted = true;
  });

  // modified from https://observablehq.com/@sw1227/reusable-2d-array-image-function
  function imshow(data, pixelSize, color) {
    // Flatten 2D input array
    const flat = [].concat.apply([], data);
    // Color Scale & Min-Max normalization
    const [min, max] = d3.extent(flat);
    const normalize = (d) => {
      const diff = max - min;
      if (Math.abs(diff) < Number.EPSILON) {
        return d;
      }
      return (d - min) / diff;
    };
    const colorScale = (d) => color(normalize(d));
    // Shape of input array
    const shape = { x: data[0].length, y: data.length };

    const context = canvas.getContext("2d");
    // Set up canvas element
    if (!canvas_initialized) {
      canvas_initialized = true;
      canvas.style.top = `${100}px`;
      canvas.style.left = `${100}px`;
      canvas.style.opacity = 1;
      canvas.style.imageRendering = "pixelated";
    }

    canvas.width = shape.x;
    canvas.height = shape.y;
    cvs_width = shape.x * pixelSize;
    cvs_height = shape.y * pixelSize;
    canvas.style.width = `${shape.x * pixelSize}px`;
    canvas.style.height = `${shape.y * pixelSize}px`;
    // const context = canvas.getContext("2d");
    // canvas.width = shape.x;
    // canvas.height = shape.y;
    // cvs_width = shape.x * pixelSize;
    // cvs_height = shape.y * pixelSize;
    // canvas.style.width = `${shape.x * pixelSize}px`;
    // canvas.style.height = `${shape.y * pixelSize}px`;
    // canvas.style.top = `${100}px`;
    // canvas.style.left = `${100}px`;
    // canvas.style.opacity = 1;
    // canvas.style.imageRendering = "pixelated";

    // Draw pixels to the canvas
    const imageData = context.createImageData(shape.x, shape.y);
    flat.forEach((d, i) => {
      let color = isNaN(d) ? { r: 0, g: 0, b: 0 } : d3.color(colorScale(d));
      imageData.data[i * 4] = color.r;
      imageData.data[i * 4 + 1] = color.g;
      imageData.data[i * 4 + 2] = color.b;
      imageData.data[i * 4 + 3] = 255;
    });
    context.putImageData(imageData, 0, 0);
  }

  let inputValue = "501";
</script>

<main>
  <canvas
    class:hide_canvas={!params.show_2d}
    bind:this={canvas}
    style="position: absolute;"
  ></canvas>
  <svg bind:this={svg} style="position: absolute; left: 0; top: 0;"></svg>
  <div id="plot" class="graph_2d"></div>
  <!-- <div id="plot_2" class="graph_2d_2"></div> -->
  <div class="checks">
    <div class="bbox" style="">
      <h4 style="display: block; margin: 0px 0px;">Idler</h4>
      <div class="controls">
        <p class="wl_label">Start</p>
        <input type="text" bind:value={params.start_idler} />
      </div>
      <div class="controls">
        <p class="wl_label">End</p>
        <input type="text" bind:value={params.end_idler} />
      </div>
    </div>

    <div class="bbox" style="">
      <h4 style="margin: 0px 0px;">Signal</h4>
      <div class="controls">
        <p class="wl_label">Start</p>
        <input type="text" bind:value={params.start_signal} />
      </div>
      <div class="controls">
        <p class="wl_label">End</p>
        <input type="text" bind:value={params.end_signal} />
      </div>
    </div>
    <div class="bbox" style="display: flex; flex-direction: column">
      <h4 style="margin: 0px 0px;">Resolution</h4>
      <div class="controls">
        <input
          type="text"
          bind:value={inputValue}
          on:input={(e) => {
            if (
              Number(e.target.value) > 100 &&
              Number(e.target.value) < 10000
            ) {
              params.numWLs = e.target.value;
            }
            if (Number(e.target.value) > 2000) {
              params.numWLs = e.target.value;
              params.show_2d = false;
            }
          }}
          class:invalid={Number(inputValue) < 100 || Number(inputValue) > 10000}
        />
      </div>
      <h4 style="margin: 0px 0px;">Show 2D</h4>
      <div class="controls">
        <label class="container">
          <input type="checkbox" bind:checked={params.show_2d} />
          <span class="checkmark"></span>
        </label>
      </div>
    </div>
  </div>
  <div class="checks">
    <label class="container"
      >Pump Envelope
      <input type="checkbox" bind:checked={params.pump_checked} />
      <span class="checkmark"></span>
    </label>

    <label class="container"
      >Phase Matching
      <input type="checkbox" bind:checked={params.phase_matching_checked} />
      <span class="checkmark"></span>
    </label>

    <label class="container"
      >Waveguide
      <input type="checkbox" bind:checked={params.waveguide_checked} />
      <span class="checkmark"></span>
    </label>
  </div>

  <div class="checks toggle">
    <p class="ppln" class:slider-off={params.ppktp_checked}>PPLN</p>
    <Toggle bind:isChecked={params.ppktp_checked} />
    <p class="ppktp" class:slider-off={!params.ppktp_checked}>PPKTP</p>
  </div>
  <div
    class="slidecontainer"
    class:slider-off={!params.pump_checked}
    style="position: relative; top: 650px;"
  >
    <p class="output">Pump WL</p>
    <input
      type="range"
      min="770"
      max="780"
      step="0.001"
      bind:value={params.WLpumpair}
      class="slider"
      id="myRange_1"
    />
    <p class="output">{params.WLpumpair} nm</p>
  </div>
  <div
    class="slidecontainer"
    class:slider-off={!params.pump_checked}
    style="position: relative; top: 650px;"
  >
    <p class="output">Pump Sigma</p>
    <input
      type="range"
      min="0.20"
      max="10"
      step="0.01"
      bind:value={params.WLSigma}
      class="slider"
      id="myRange_1"
    />
    <p class="output">{params.WLSigma} nm</p>
  </div>
  <div
    class="slidecontainer"
    class:slider-off={!params.waveguide_checked}
    style="position: relative; top: 650px;"
  >
    <p class="output">Cavity Length</p>
    <input
      type="range"
      min=".999"
      max="1.001"
      step="0.000001"
      bind:value={params.L}
      class="slider"
      id="myRange_1"
    />
    <p
      style="font-size: 1.8rem; padding-left: 10px; padding-right: 10px; margin:0; margin-top: 7px;"
    >
      ×
    </p>
    <input
      type="number"
      style="height: 29px; width: 70px; margin-top: 12px;"
      bind:value={params.L_multiplier}
    />
    <p class="output">
      {(params.L * params.L_multiplier).toFixed(4)} mm
    </p>
    <input
      type="number"
      style="height: 29px; width: 70px; margin-top: 12px;"
      bind:value={params.L}
    />
  </div>
  <div
    class="slidecontainer"
    class:slider-off={!(params.waveguide_checked || params.phase_matching_checked)}
    style="position: relative; top: 650px;"
  >
    <p class="output">Temperature</p>
    <input
      type="range"
      min="-50"
      max="100"
      step="1"
      bind:value={params.Texpt}
      class="slider"
      id="myRange_1"
    />
    <p class="output">{params.Texpt} Celsius</p>
  </div>
  <div
    class="slidecontainer"
    class:slider-off={!(params.waveguide_checked || params.phase_matching_checked)}
    style="position: relative; top: 650px;"
  >
    <p class="output">Polling Period</p>
    <input
      type="range"
      min="20"
      max="30"
      step="0.001"
      bind:value={params.polling_period}
      class="slider"
      id="myRange_1"
    />
    <p class="output">{params.polling_period} um</p>
  </div>
  <div
    class="slidecontainer"
    class:slider-off={!params.waveguide_checked}
    style="position: relative; top: 650px;"
  >
    <p class="output">R1</p>
    <input
      type="range"
      min="0"
      max="1"
      step=".0001"
      bind:value={params.R1}
      class="slider"
      id="myRange_1"
    />
    <p class="output">{params.R1}</p>
  </div>
  <div
    class="slidecontainer"
    class:slider-off={!params.waveguide_checked}
    style="position: relative; top: 650px;"
  >
    <p class="output">R2</p>
    <input
      type="range"
      min="0"
      max="1"
      step=".0001"
      bind:value={params.R2}
      class="slider"
      id="myRange_1"
    />
    <p class="output">{params.R2}</p>
  </div>
  <div class="footer"></div>
</main>

<style>
  .hide_canvas {
    opacity: 0.1;
  }

  .invalid {
    border: 1px solid red;
  }

  .wl_label {
    font-size: 1.2rem;
    padding: 5px 5px;
    margin: 0;
    margin-top: 7px;
  }

  .controls {
    display: flex;
    flex-direction: row;
    justify-content: space-between;
    margin: 10px 10px;
    /* padding-left: 10px;
    padding-right: 10px; */
  }

  .footer {
    height: 800px;
  }
  .graph_2d {
    position: relative;
    top: 90px;
    left: 740px;
    width: 1500px;
    height: 0px;
  }

  .graph_2d_2 {
    position: relative;
    top: 90px;
    left: 740px;
    width: 1500px;
    height: 1500px;
  }

  .ppln,
  .ppktp {
    font-size: 1.3rem;
    padding-left: 10px;
    padding-right: 10px;
    margin: 0;
    margin-top: 7px;
  }

  .slider-off {
    opacity: 0.2;
  }

  .deactivate {
    color: rgb(215, 215, 215);
  }
  .checks {
    position: relative;
    top: 650px;
    width: 700px; /* Width of the outside container */
    display: flex;
    flex-direction: row;
    justify-content: space-around;
    /* padding-right: 50px; */
    box-shadow: 0 0 12px rgba(0, 0, 0, 0.07);
    border: 1.3px solid #e5e5e5;
    padding-top: 0.8rem;
    margin-bottom: 1rem;
  }

  .toggle {
    padding-bottom: 0.8rem;
    justify-content: space-evenly;
  }

  /* Customize the label (the container) */
  .container {
    display: block;
    position: relative;
    padding-left: 35px;
    margin-bottom: 12px;
    cursor: pointer;
    font-size: 22px;
    -webkit-user-select: none;
    -moz-user-select: none;
    -ms-user-select: none;
    user-select: none;
  }

  /* Hide the browser's default checkbox */
  .container input {
    position: absolute;
    opacity: 0;
    cursor: pointer;
    height: 0;
    width: 0;
  }

  /* Create a custom checkbox */
  .checkmark {
    position: absolute;
    top: 0;
    left: 0;
    height: 25px;
    width: 25px;
    background-color: #eee;
  }

  /* On mouse-over, add a grey background color */
  .container:hover input ~ .checkmark {
    background-color: #ccc;
  }

  /* When the checkbox is checked, add a blue background */
  .container input:checked ~ .checkmark {
    background-color: #2196f3;
  }

  /* Create the checkmark/indicator (hidden when not checked) */
  .checkmark:after {
    content: "";
    position: absolute;
    display: none;
  }

  /* Show the checkmark when checked */
  .container input:checked ~ .checkmark:after {
    display: block;
  }

  /* Style the checkmark/indicator */
  .container .checkmark:after {
    left: 9px;
    top: 5px;
    width: 5px;
    height: 10px;
    border: solid white;
    border-width: 0 3px 3px 0;
    -webkit-transform: rotate(45deg);
    -ms-transform: rotate(45deg);
    transform: rotate(45deg);
  }

  .slidecontainer {
    width: 700px; /* Width of the outside container */
    /* padding-top: 1rem; */
    /* background-color: #04aa6d; */
    display: flex;
    flex-direction: row;
    margin-left: 10px;
  }

  /* The slider itself */
  .slider {
    -webkit-appearance: none; /* Override default CSS styles */
    appearance: none;
    width: 50%; /* Full-width */
    height: 25px;
    background: #e5e5e5; /* Grey background */
    outline: none; /* Remove outline */
    /* opacity: 0.7;  */
    /* -webkit-transition: 0.2s;  */
    /* transition: opacity 0.2s; */
    border-radius: 0.3rem;
    margin-top: 13px;
  }

  /* Mouse-over effects */
  /* .slider:hover {
    opacity: 1; 
  }  */

  /* The slider handle (use -webkit- (Chrome, Opera, Safari, Edge) and -moz- (Firefox) to override default look) */
  .slider::-webkit-slider-thumb {
    -webkit-appearance: none; /* Override default look */
    appearance: none;
    width: 25px; /* Set a specific slider handle width */
    height: 25px; /* Slider handle height */
    background: #9a9a9a; /* Green background */
    cursor: pointer; /* Cursor on hover */
    border-radius: 0.3rem;
  }

  .slider::-moz-range-thumb {
    width: 25px; /* Set a specific slider handle width */
    height: 25px; /* Slider handle height */
    background: #04aa6d; /* Green background */
    cursor: pointer; /* Cursor on hover */
  }

  .output {
    /* font style */
    font-family: "Roboto", sans-serif;
    font-weight: bold;
    /* padding-bottom: 10px; */
    padding-left: 10px;
    padding-right: 10px;
    min-width: 140px;
  }

  :global(svg g path),
  :global(svg g line) {
    fill: none;
    stroke: black;
    stroke-width: 1.3px; /* Make the axis lines thicker */
  }

  :global(svg g text) {
    font-size: 15px; /* Make the text bigger */
  }
</style>
