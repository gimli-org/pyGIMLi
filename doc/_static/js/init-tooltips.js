document.addEventListener("DOMContentLoaded", function () {
  if (typeof bootstrap === "undefined" || !bootstrap.Tooltip) return;
  document
    .querySelectorAll('[data-bs-toggle="tooltip"]')
    .forEach(function (el) {
      new bootstrap.Tooltip(el);
    });
});
