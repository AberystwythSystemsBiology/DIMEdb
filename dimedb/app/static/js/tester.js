var app = angular.module("myApp", []);


function renderMetabolite(metabolite, $scope) {
  $scope.metaboliteName = metabolite["Identification Information"]["Name"];
}

app.controller("MetaboliteViewer", function($scope, $http) {
  $scope.loadPeople = function() {
    $http.get("http://localhost:5000/api/metabolites?where={%22_id%22%20:%20%22WFWLQNSHRPWKFK-UHFFFAOYSA-N%22}")
      .then(function(data) {
        renderMetabolite(data["data"]["_items"][0], $scope);
      })
  };
  $scope.loadPeople();
});
