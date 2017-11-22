webApp = angular.module("DIMEdb", []);

webApp.controller("NavbarController", ["$scope", function($scope) {
  $(".navbar-brand img").hide();
}]);
