webApp = angular.module("DIMEdb", []);

webApp.controller("NavbarController", ["$scope", function($scope) {}]);

// Retrieves metabolite information from API.
webApp.factory("getMetabolite", function($http) {
  data = $http.get(window.location + "/json").then(function(response) {
    return response.data;
  });
  return data
});

function getMolecularMatches(molecularFormula) {
  var api_url = getBaseURL() + "api/metabolites/?where={";
  api_url += '"Identification Information.Molecular Formula" : "';
  api_url += molecularFormula + '"}';
  api_url += '&projection={"Identification Information.Name" : 1}';

  var metabolites = get_metabolites(api_url);

  $("#molecularMatchTable > tbody").empty();

  for(idx in metabolites) {
    var metabolite = metabolites[idx];
    var name = metabolite["Identification Information"]["Name"];
    var inchikey = metabolite["_id"];

    var structure_url = getBaseURL() + "view/" + inchikey;
    var image_url = structure_url+"/structure";

    $("#molecularMatchTable > tbody").append(
      "<tr><td><img class='img-responsive img-circle'src='"+image_url+"'></td><td>"
              +"<a href='"+structure_url+"' target='_blank'>"+name+"</a></td></tr>");
  }

  $('#molecularMatchTable').dataTable( {
    "columnDefs": [
      { "width": "10%", "targets": 0 }
    ]
  });
}

// Take in the adduct info and polarity, fill the selector with info.
function generateAdductSelection(adduct_info, polarity) {
  $("#adduct_selector").empty();
  for (var idx in adduct_info) {
    var adduct = adduct_info[idx];
    if (adduct["Polarity"] == polarity) {
      $("#adduct_selector").append($("<option>", {
        value: adduct["Adduct"],
        text: adduct["Adduct"] + ": " + adduct["Accurate Mass"] + " m/z"
      }));
    }
  }

  if ($("#adduct_selector").length < 1) {
    $("#adduct_selector").prop("disabled", true);
  } else {
    $("#adduct_selector").prop("disabled", false);
  }
}

function getMassesAndIntensities(adduct_info, polarity, adduct_val) {

  var masses = [];
  var intensities = [];

  for (var idx in adduct_info) {
    var adduct = adduct_info[idx]
    if (adduct["Polarity"] == polarity && adduct["Adduct"] == adduct_val) {
      var isotopic_distribution = adduct["Isotopic Distribution"];
      for (d_idx in isotopic_distribution) {
        d = isotopic_distribution[d_idx]
        masses.push(d[0]);
        intensities.push(d[1]);
        }
      }
    }
  return [masses, intensities]
}

function generateChartAndTable(masses, intensities, chartHeight, chartWidth) {
  // If this is true, limit values to > 5%.
  var limited = $("#rel_i_limiter").prop('checked');

  // Table Filler
  $("#distribution_table tbody > tr").remove();

  for(idx in masses) {
    var mass = masses[idx];
    var intensity = intensities[idx];
    if (intensity >= 5) {
      $("#distribution_table > tbody").append("<tr><td><b>"+mass+"</b></td><td><b>"
                +intensity+"</b></td></tr>");
    }
    else {
      $("#distribution_table > tbody").append("<tr><td>"+mass+"</td><td>"
                +intensity+"</td></tr>");
    }
  }


  if (limited == true) {
    var m = [];
    var i = [];
    for (idx in masses) {
      var mass = masses[idx];
      var intensity = intensities[idx];
      if (intensity >= 5) {
        m.push(mass);
        i.push(intensity);
      }
    }

    var masses = m;
    var intensities = i;
  }

  function plot() {
    var d3 = Plotly.d3;

    var gd3 = d3.select("#myChart");

    var gd = gd3.node();

    var trace1 = {
      x: masses,
      y: intensities,
      marker: {color: 'rgb(55, 83, 109)'},
      type: 'bar'
    };

    var data = [trace1];

    if (masses.length > 3) {
      var bg = 0.90
    }
    else {
      var bg = 0.99
    }

    var layout={
      height: chartHeight,
      width: chartWidth,
      margin: {
      b: 35,
      t: 10
    },
      bargap: bg,
      xaxis : {
        "title" : "mass-to-ion (m/z)"
      },
      yaxis : {
        "title" : "relative abundance (%)"
      },
      showlegend: false,
      autosize: true
  }

    Plotly.newPlot(gd, data, layout, {displayModeBar: false});

  }

  plot();

}

function generateChart($scope, refresh_adduct=true) {
  var polarity = $("input[name=ionisation]:checked").val();
  if (refresh_adduct == true) {
    generateAdductSelection($scope.metabolite["Adducts"], polarity);
  }
  var adduct_val = $("#adduct_selector").val();
  var values = getMassesAndIntensities($scope.metabolite["Adducts"], polarity, adduct_val);

  var chartHeight = $("#myChart").height();

  var chartWidth = $("#myChart").width();

  generateChartAndTable(values[0], values[1], chartHeight, chartWidth);
}


function getSkeletons(inchikey) {
    var api_url = getBaseURL() + "api/metabolites/?where={";
    api_url += '"_id" : { "$regex" : "^' + inchikey.split("-")[0] +'"}}';
    api_url += '&projection={"Identification Information.Name" : 1}';

    var metabolites = get_metabolites(api_url);

    var any = false;

    for (index in metabolites) {
        var m = metabolites[index];
        if (m["_id"] != inchikey) {
           $("#skeleton_list").append("<a href='" + getBaseURL() +"view/" + m["_id"] + "' target='_blank'><li class='list-group-item'>"
            +m["Identification Information"]["Name"] + "</li></a>");
            any = true;
        }
    }

    if (any == false) {
        $("#skeleton_list").append("<li class='list-group-item'> None available </li>");
    }
}

function fillJCAMPTool($scope) {
  var negative = [];
  var positive = [];
  var neutral = [];
  for(idx in $scope.metabolite["Adducts"]) {
    var adduct_info = $scope.metabolite["Adducts"][idx];
    if(adduct_info["Polarity"] == "Negative") {
      negative.push(adduct_info)
    }
    else if (adduct_info["Polarity"] == "Positive") {
      positive.push(adduct_info);
    }
    else {
      neutral.push(adduct_info);
    }
  }

  console.log(negative, positive, neutral);
}

webApp.controller("MetaboliteView", function(getMetabolite, $scope) {
  // Generate inchikey from URL
  var inchikey = window.location.pathname.split("/")[2];

  getMetabolite.then(function(success) {
    $scope.metabolite = success

    // Change molecular formula to standard form
    var molecularFormula = $("#molecularFormula").html();
    $("#molecularFormula").html(molecularFormula.replace(/([0-9]+)/g, '<sub>$1</sub>'));

    var inchikey = $("#inchikey").html();
    getSkeletons(inchikey);
    // Generate chart for Neutral
    generateChart($scope);


    $("#loadingPanel").fadeOut(100);

    $("#metaboliteDiv").fadeIn(500);

    $("input[type=radio][name=ionisation]").change(function() {
      generateChart($scope);
    });

    $("#adduct_selector").change(function() {
      generateChart($scope, false);
    });

    $("#rel_i_limiter").change(function() {
      generateChart($scope, false);
    })

    $(window).resize(function() {
      generateChart($scope, false);
    });

    $("#formulaSearchBtn").click(function() {
      getMolecularMatches(molecularFormula);
      $("#formula_search_modal").modal("toggle");
    });
  });

});
