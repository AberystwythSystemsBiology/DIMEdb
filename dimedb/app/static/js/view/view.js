webApp = angular.module("DIMEdb", []);

webApp.controller("NavbarController", ["$scope", function($scope) {}]);

// Retrieves metabolite information from API.
webApp.factory("getMetabolite", function($http) {
  data = $http.get(window.location + "/json").then(function(response) {
    return response.data;
  });
  return data
});

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
  var limited = $("#rel_i_limiter").prop('checked');

  var masses = [];
  var intensities = [];

  for (var idx in adduct_info) {
    var adduct = adduct_info[idx]
    if (adduct["Polarity"] == polarity && adduct["Adduct"] == adduct_val) {
      var isotopic_distribution = adduct["Isotopic Distribution"];
      for (d_idx in isotopic_distribution) {
        d = isotopic_distribution[d_idx]
        if (limited == false) {
          masses.push(d[0]);
          intensities.push(d[1]);
        } else {
          if (d[1] > 5.0) {
            masses.push(d[0]);
            intensities.push(d[1]);
          }
        }
      }
    }
  }
  return [masses, intensities]
}

function generateChartAndTable(masses, intensities) {
  var chart = $("#myChart");

  var myChart = new Chart(chart, {
    type: 'line',
    data: {
      labels: masses,
      datasets: [{
        data: intensities
      }]
    },
    options: {
      responsive: true,
      maintainAspectRatio: false,
      legend: {
        display: false
      },
      scales: {
        xAxes: [{
          display: true,
          ticks : {

          },
          scaleLabel: {
            display: true,
            labelString: "Mass-to-ion (m/z)"
          }
        }],
        yAxes: [{
          display: true,
          ticks : {
            beginAtZero: true,
            stepSize : 25
          },
          scaleLabel: {
            display: true,
            labelString: "Relative Intensity (%)"
          }
        }]
      }
    }
  });
}

function generateChart($scope, refresh_adduct=true) {
  var polarity = $("input[name=ionisation]:checked").val();
  if (refresh_adduct == true) {
    generateAdductSelection($scope.metabolite["Adducts"], polarity);
  }
  var adduct_val = $("#adduct_selector").val();
  var values = getMassesAndIntensities($scope.metabolite["Adducts"], polarity, adduct_val);
  generateChartAndTable(values[0], values[1]);
}

webApp.controller("MetaboliteView", function(getMetabolite, $scope) {
  // Generate inchikey from URL
  var inchikey = window.location.pathname.split("/")[2];

  getMetabolite.then(function(success) {
    $scope.metabolite = success

    // Change molecular formula to standard form
    var molecularFormula = $("#molecularFormula").html();
    $("#molecularFormula").html(molecularFormula.replace(/([0-9]+)/g, '<sub>$1</sub>'));

    // Generate chart for Neutral
    generateChart($scope);

    $("input[type=radio][name=ionisation]").change(function() {
      generateChart($scope);
    });

    $("#adduct_selector").change(function() {
      generateChart($scope, false);
    });

    $("#rel_i_limiter").change(function() {
      generateChart($scope, false);
    })
  });

});
