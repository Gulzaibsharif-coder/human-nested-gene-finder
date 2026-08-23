const state = {
  allData: [],
  filtered: [],
  page: 1,
  pageSize: 25,
  currentRow: null
};

const els = {
  searchInput: document.getElementById("searchInput"),
  searchButton: document.getElementById("searchButton"),
  resetButton: document.getElementById("resetButton"),
  chromosomeFilter: document.getElementById("chromosomeFilter"),
  hostTypeFilter: document.getElementById("hostTypeFilter"),
  nestedTypeFilter: document.getElementById("nestedTypeFilter"),
  directionFilter: document.getElementById("directionFilter"),
  resultCount: document.getElementById("resultCount"),
  resultsBody: document.getElementById("resultsBody"),
  prevPage: document.getElementById("prevPage"),
  nextPage: document.getElementById("nextPage"),
  pageInfo: document.getElementById("pageInfo"),
  downloadButton: document.getElementById("downloadButton"),
  detailsModal: document.getElementById("detailsModal"),
  modalBackdrop: document.getElementById("modalBackdrop"),
  modalClose: document.getElementById("modalClose"),
  modalTitle: document.getElementById("modalTitle"),
  modalContent: document.getElementById("modalContent"),
  copyIdButton: document.getElementById("copyIdButton"),
  downloadRowButton: document.getElementById("downloadRowButton")
};

function escapeHtml(value) {
  return String(value ?? "")
    .replaceAll("&", "&amp;")
    .replaceAll("<", "&lt;")
    .replaceAll(">", "&gt;")
    .replaceAll('"', "&quot;")
    .replaceAll("'", "&#039;");
}

function formatNumber(value) {
  return Number(value).toLocaleString();
}

function coordinateText(start, end) {
  return `${formatNumber(start)}–${formatNumber(end)}`;
}

function uniqueSorted(values) {
  return [...new Set(values.filter(Boolean))]
    .sort((a, b) => String(a).localeCompare(String(b), undefined, { numeric: true }));
}

async function loadData() {
  const [dataResponse, statsResponse] = await Promise.all([
    fetch("data/nested_genes.json"),
    fetch("data/stats.json")
  ]);

  if (!dataResponse.ok || !statsResponse.ok) {
    throw new Error("Unable to load dataset files.");
  }

  state.allData = await dataResponse.json();
  const stats = await statsResponse.json();

  document.getElementById("heroRelationships").textContent = formatNumber(stats.relationships);
  document.getElementById("heroHosts").textContent = formatNumber(stats.host_genes);
  document.getElementById("heroNested").textContent = formatNumber(stats.nested_genes);

  document.getElementById("statGenes").textContent = formatNumber(stats.genes);
  document.getElementById("statTranscripts").textContent = formatNumber(stats.transcripts);
  document.getElementById("statExons").textContent = formatNumber(stats.exons);
  document.getElementById("statIntrons").textContent = formatNumber(stats.introns);

  populateFilters();
  applyFilters();
}

function populateFilters() {
  uniqueSorted(state.allData.map(row => row.chromosome))
    .forEach(value => {
      const option = document.createElement("option");
      option.value = value;
      option.textContent = value;
      els.chromosomeFilter.appendChild(option);
    });

  uniqueSorted(state.allData.map(row => row.host_type))
    .forEach(value => {
      const option = document.createElement("option");
      option.value = value;
      option.textContent = value;
      els.hostTypeFilter.appendChild(option);
    });

  uniqueSorted(state.allData.map(row => row.nested_type))
    .forEach(value => {
      const option = document.createElement("option");
      option.value = value;
      option.textContent = value;
      els.nestedTypeFilter.appendChild(option);
    });
}

function parseCoordinateQuery(query) {
  const match = query.match(/^([^:]+):(\d+)-(\d+)$/);
  if (!match) return null;

  return {
    chromosome: match[1],
    start: Number(match[2]),
    end: Number(match[3])
  };
}

function rowMatchesCoordinate(row, coordinate) {
  if (!coordinate) return false;
  if (row.chromosome !== coordinate.chromosome) return false;

  return (
    Number(row.nested_start) <= coordinate.end &&
    Number(row.nested_end) >= coordinate.start
  );
}

function matchesSearch(row, query) {
  if (!query) return true;

  const coordinate = parseCoordinateQuery(query);

  if (coordinate) {
    return rowMatchesCoordinate(row, coordinate);
  }

  const q = query.toLowerCase();

  const fields = [
    row.host_id,
    row.host_name,
    row.host_type,
    row.nested_id,
    row.nested_name,
    row.nested_type,
    row.chromosome
  ];

  return fields.some(value =>
    String(value ?? "").toLowerCase().includes(q)
  );
}

function applyFilters() {
  const query = els.searchInput.value.trim();
  const chromosome = els.chromosomeFilter.value;
  const hostType = els.hostTypeFilter.value;
  const nestedType = els.nestedTypeFilter.value;
  const direction = els.directionFilter.value;

  state.filtered = state.allData.filter(row => {

    if (!matchesSearch(row, query)) return false;

    if (chromosome && row.chromosome !== chromosome) return false;

    if (hostType && row.host_type !== hostType) return false;

    if (nestedType && row.nested_type !== nestedType) return false;

    if (direction === "host" && !row.host_id) return false;
    if (direction === "nested" && !row.nested_id) return false;

    return true;
  });

  state.page = 1;
  renderResults();
}

function renderResults() {
  const total = state.filtered.length;
  const totalPages = Math.max(1, Math.ceil(total / state.pageSize));

  if (state.page > totalPages) {
    state.page = totalPages;
  }

  const start = (state.page - 1) * state.pageSize;
  const rows = state.filtered.slice(start, start + state.pageSize);

  els.resultCount.textContent =
    `${formatNumber(total)} matching relationship${total === 1 ? "" : "s"}`;

  els.pageInfo.textContent = `Page ${state.page} of ${totalPages}`;

  els.prevPage.disabled = state.page <= 1;
  els.nextPage.disabled = state.page >= totalPages;

  if (rows.length === 0) {
    els.resultsBody.innerHTML = `
      <tr>
        <td colspan="8" class="empty-cell">
          No matching relationships found.
        </td>
      </tr>
    `;
    return;
  }

  els.resultsBody.innerHTML = rows.map((row, index) => `
    <tr>
      <td>
        <strong>${escapeHtml(row.nested_name || row.nested_id)}</strong><br>
        <span>${escapeHtml(row.nested_id)}</span>
      </td>

      <td>
        <span class="type-pill">${escapeHtml(row.nested_type)}</span>
      </td>

      <td>
        <strong>${escapeHtml(row.host_name || row.host_id)}</strong><br>
        <span>${escapeHtml(row.host_id)}</span>
      </td>

      <td>
        <span class="type-pill">${escapeHtml(row.host_type)}</span>
      </td>

      <td>${escapeHtml(row.chromosome)}</td>

      <td>
        ${coordinateText(row.nested_start, row.nested_end)}
      </td>

      <td>
        ${coordinateText(row.intron_start, row.intron_end)}
      </td>

      <td>
        <button class="table-action" data-index="${index}">
          Details
        </button>
      </td>
    </tr>
  `).join("");

  [...els.resultsBody.querySelectorAll(".table-action")]
    .forEach(button => {
      button.addEventListener("click", () => {
        const row = rows[Number(button.dataset.index)];
        openModal(row);
      });
    });
}

function openModal(row) {
  state.currentRow = row;

  els.modalTitle.textContent =
    row.nested_name || row.nested_id;

  els.modalContent.innerHTML = `
    ${detail("Nested gene", row.nested_name || row.nested_id)}
    ${detail("Nested Ensembl ID", row.nested_id)}
    ${detail("Nested gene type", row.nested_type)}
    ${detail("Host gene", row.host_name || row.host_id)}
    ${detail("Host Ensembl ID", row.host_id)}
    ${detail("Host gene type", row.host_type)}
    ${detail("Chromosome", row.chromosome)}
    ${detail("Nested coordinates", coordinateText(row.nested_start, row.nested_end))}
    ${detail("Containing intron", coordinateText(row.intron_start, row.intron_end))}
    ${detail("Containment", "Fully contained")}
  `;

  els.detailsModal.classList.remove("hidden");
}

function detail(label, value) {
  return `
    <div class="detail-item">
      <span class="detail-label">${escapeHtml(label)}</span>
      <span class="detail-value">${escapeHtml(value)}</span>
    </div>
  `;
}

function closeModal() {
  els.detailsModal.classList.add("hidden");
  state.currentRow = null;
}

function downloadCsv(rows, filename) {
  const header = [
    "host_id",
    "host_name",
    "host_type",
    "nested_id",
    "nested_name",
    "nested_type",
    "chromosome",
    "nested_start",
    "nested_end",
    "intron_start",
    "intron_end"
  ];

  const csvRows = [
    header,
    ...rows.map(row =>
      header.map(field => {
        const value = row[field] ?? "";
        return `"${String(value).replaceAll('"', '""')}"`;
      })
    )
  ];

  const blob = new Blob(
    [csvRows.map(row => row.join(",")).join("\n")],
    { type: "text/csv;charset=utf-8;" }
  );

  const url = URL.createObjectURL(blob);
  const link = document.createElement("a");
  link.href = url;
  link.download = filename;
  link.click();
  URL.revokeObjectURL(url);
}

els.searchButton.addEventListener("click", applyFilters);

els.searchInput.addEventListener("keydown", event => {
  if (event.key === "Enter") {
    applyFilters();
  }
});

[
  els.chromosomeFilter,
  els.hostTypeFilter,
  els.nestedTypeFilter,
  els.directionFilter
].forEach(control => control.addEventListener("change", applyFilters));

els.resetButton.addEventListener("click", () => {
  els.searchInput.value = "";
  els.chromosomeFilter.value = "";
  els.hostTypeFilter.value = "";
  els.nestedTypeFilter.value = "";
  els.directionFilter.value = "both";
  applyFilters();
});

document.querySelectorAll(".example-query").forEach(button => {
  button.addEventListener("click", () => {
    els.searchInput.value = button.dataset.query;
    applyFilters();
    document.getElementById("results").scrollIntoView({ behavior: "smooth" });
  });
});

els.prevPage.addEventListener("click", () => {
  if (state.page > 1) {
    state.page--;
    renderResults();
  }
});

els.nextPage.addEventListener("click", () => {
  const totalPages = Math.ceil(state.filtered.length / state.pageSize);

  if (state.page < totalPages) {
    state.page++;
    renderResults();
  }
});

els.downloadButton.addEventListener("click", () => {
  downloadCsv(state.filtered, "nested_gene_results.csv");
});

els.modalClose.addEventListener("click", closeModal);
els.modalBackdrop.addEventListener("click", closeModal);

document.addEventListener("keydown", event => {
  if (event.key === "Escape") {
    closeModal();
  }
});

els.copyIdButton.addEventListener("click", async () => {
  if (!state.currentRow) return;

  try {
    await navigator.clipboard.writeText(state.currentRow.nested_id);
    els.copyIdButton.textContent = "Copied ✓";

    setTimeout(() => {
      els.copyIdButton.textContent = "Copy nested ID";
    }, 1400);
  } catch {
    els.copyIdButton.textContent = "Copy unavailable";
  }
});

els.downloadRowButton.addEventListener("click", () => {
  if (!state.currentRow) return;

  downloadCsv(
    [state.currentRow],
    `${state.currentRow.nested_id}_nested_relationship.csv`
  );
});

loadData().catch(error => {
  console.error(error);

  els.resultCount.textContent =
    "Unable to load the dataset. Run the site from a local web server or GitHub Pages.";

  els.resultsBody.innerHTML = `
    <tr>
      <td colspan="8" class="empty-cell">
        Dataset could not be loaded.
      </td>
    </tr>
  `;
});
